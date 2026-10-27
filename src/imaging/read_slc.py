#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""read_slc.py — Sentinel-1 Level-1 SLC 读取器（A 组 · 第2周另一条线）

读取 SLC 测量 TIFF 的复数图像，解析 annotation XML 成像参数，可选 ESA 辐射定标，
抽取 ROI 训练样本，保存复数 .npy + 参数 .json（+ 可选 MATLAB .mat + 幅度预览 PNG）。

本产品 SLC 测量 TIFF 实测格式：
  - SampleFormat=5（complex）、SamplesPerPixel=1、BitsPerSample=32
  - 实际为 **int16 复数（I/Q 交织，每像素 4 字节）**、小端、无压缩
  - 分块方式：本产品 RowsPerStrip=1（逐行一个 strip）
  整景为 39297×20104 ≈ 6.3 GB complex64，直接整景载入会撑爆内存；本模块按
  “行 strip”**按需读取所需块**（只读目标行、并在行内再切列），因此可低成本取任意块。

依赖：numpy, tifffile（已装入 data/raw/Sentinel-L0条带/.venv）；
      可选：scipy（写 .mat / 定标插值）、Pillow（幅度预览 PNG）。

典型命令行：
  python read_slc.py --safe_dir "<SLC.SAFE>" --polarization VV
作为库：
  from read_slc import read_slc
  slc, params = read_slc(SAFE, polarization='VV', row0=6000, row1=8048, col0=6000, col1=8048)
"""
from __future__ import annotations

import os
import glob
import json
import struct
import argparse
from datetime import datetime

import numpy as np

C = 299792458.0  # 光速 m/s


# --------------------------------------------------------------------------
# 1. 文件定位
# --------------------------------------------------------------------------
def find_pol_files(safe_dir, polarization):
    """返回 (annotation_xml, measurement_tiff, calibration_xml|None)。"""
    pol = polarization.lower()
    ann = glob.glob(os.path.join(safe_dir, "annotation", f"*-slc-{pol}-*.xml"))
    meas = glob.glob(os.path.join(safe_dir, "measurement", f"*-slc-{pol}-*.tiff"))
    cal = glob.glob(os.path.join(safe_dir, "annotation", "calibration", f"*-slc-{pol}-*.xml"))
    assert ann, f"未找到极化 {pol} 的 annotation XML"
    assert meas, f"未找到极化 {pol} 的测量 TIFF"
    return ann[0], meas[0], (cal[0] if cal else None)


# --------------------------------------------------------------------------
# 2. annotation XML 参数解析
# --------------------------------------------------------------------------
def _first_text(root, tag):
    els = list(root.iter(tag))
    return (els[0].text or "").strip() if els else None


def _float(txt):
    try:
        return float(txt)
    except (TypeError, ValueError):
        return None


def parse_annotation(xml_path):
    import xml.etree.ElementTree as ET

    r = ET.parse(xml_path).getroot()
    p = {}

    str_tags = ["missionId", "productType", "mode", "swath", "polarisation",
                "startTimeUTC", "stopTimeUTC",
                "productFirstLineUtcTime", "productLastLineUtcTime", "pass"]
    int_tags = ["numberOfSamples", "numberOfLines", "absoluteOrbitNumber", "missionDataTakeId"]
    flt_tags = ["slantRangeTime", "rangePixelSpacing", "azimuthPixelSpacing",
                "azimuthTimeInterval", "azimuthFrequency", "incidenceAngleMidSwath",
                "zeroDopMinusAcqTime", "platformHeading", "rangeSamplingRate",
                "radarFrequency", "azimuthSteeringRate"]
    for tag in str_tags:
        v = _first_text(r, tag)
        if v is not None:
            p[tag] = v
    for tag in int_tags:
        v = _first_text(r, tag)
        if v is not None:
            p[tag] = int(float(v))
    for tag in flt_tags:
        v = _first_text(r, tag)
        if v is not None:
            p[tag] = _float(v)

    # 通用补充检索（标签可能随产品模板略不同）
    for tag in ["rangeBandwidth", "azimuthBandwidth", "rangeResolution",
                "azimuthResolution", "sceneCenterLatitude", "sceneCenterLongitude"]:
        if tag not in p:
            v = _first_text(r, tag)
            if v is not None:
                p[tag] = _float(v)

    p["orbit"] = _parse_orbit(r)

    # 派生量
    fc = p.get("radarFrequency")
    if fc:
        p["carrier_frequency_hz"] = fc
        p["wavelength_m"] = C / fc
    if "azimuthFrequency" in p:
        p["prf_hz"] = p["azimuthFrequency"]
    return p


def _parse_orbit(root):
    """提取 orbitList 全部状态向量（UTC 时间 + ECEF 位置/速度）。"""
    orbits = list(root.iter("orbit"))
    if not orbits:
        return {"available": False, "reason": "annotation 中无 <orbit>"}
    times, pos, vel = [], [], []
    for o in orbits:
        times.append((o.findtext("time") or "").strip())

        def xyz(e):
            if e is None:
                return None
            return [_float(e.findtext("x")), _float(e.findtext("y")), _float(e.findtext("z"))]

        pos.append(xyz(o.find("position")))
        vel.append(xyz(o.find("velocity")))
    return {"available": True, "count": len(times),
            "times_utc": times,
            "position_ecef_m": pos,
            "velocity_ecef_mps": vel}


# --------------------------------------------------------------------------
# 3. SLC TIFF 只读 / 分块读取
# --------------------------------------------------------------------------
class SlcTiff:
    """解析 SLC 测量 TIFF 的 IFD，并按行 strip 按需读取复数块（避免整景载入）。"""

    def __init__(self, tiff_path):
        self.path = tiff_path
        self._parse_ifd()

    def _parse_ifd(self):
        with open(self.path, "rb") as f:
            hdr = f.read(8)
            self.endian = "<" if hdr[:2] == b"II" else ">"
            ifd_off = struct.unpack(self.endian + "I", hdr[4:8])[0]
            f.seek(ifd_off)
            n = struct.unpack(self.endian + "H", f.read(2))[0]
            raw = f.read(n * 12)

        SZX = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2, 9: 4, 10: 8, 11: 4, 12: 8}
        ent = {}
        for i in range(n):
            o = i * 12
            tag, typ, cnt = struct.unpack(self.endian + "HHI", raw[o:o + 8])
            ent[tag] = (typ, cnt, raw[o + 8:o + 12])

        def val(t):
            if t not in ent:
                return None
            typ, cnt, vo = ent[t]
            tot = SZX[typ] * cnt
            if tot <= 4:
                b = vo[:tot]
            else:
                off = struct.unpack(self.endian + "I", vo[:4])[0]
                with open(self.path, "rb") as g:
                    g.seek(off)
                    b = g.read(tot)
            if typ in (3, 8):
                return list(struct.unpack(self.endian + "%dH" % cnt, b))
            if typ in (4, 9):
                return list(struct.unpack(self.endian + "%dI" % cnt, b))
            return b

        self.width = val(256)[0]
        self.length = val(257)[0]
        bits = val(258)
        self.bits = bits[0] if isinstance(bits, list) else None
        self.sample_format = val(339)[0] if val(339) else 1
        self.samples_per_pixel = val(277)[0]
        self.compression = val(259)[0]
        self.planar = val(284)[0] if val(284) else 1
        self.rows_per_strip = val(278)[0] if val(278) else self.length
        self.strip_offsets = np.array(val(273), dtype=np.int64)
        self.strip_byte_counts = np.array(val(279), dtype=np.int64)

        bpp = int(self.strip_byte_counts[0] // (self.rows_per_strip * self.width))
        self.bytes_per_pixel = bpp
        if bpp == 4:
            self.sample_dtype = np.int16
        elif bpp == 8:
            self.sample_dtype = np.float32
        else:
            raise ValueError(f"不支持的 SLC 像元字节数 {bpp}（仅支持 4/8）")

        if self.compression != 1:
            raise ValueError(
                f"仅支持无压缩(Compression=1) SLC，当前 Compression={self.compression}；"
                f"可用 tifffile 整景读取作为替代。")
        if self.sample_format != 5:
            raise ValueError(f"仅支持 complex(SampleFormat=5)，当前 {self.sample_format}")

    def read_block(self, row0=0, row1=None, col0=0, col1=None):
        if row1 is None:
            row1 = self.length
        if col1 is None:
            col1 = self.width
        row0 = max(0, int(row0))
        row1 = min(self.length, int(row1))
        col0 = max(0, int(col0))
        col1 = min(self.width, int(col1))
        if row1 <= row0 or col1 <= col0:
            raise ValueError("读取窗口无效：row1>row0 且 col1>col0")

        rps = self.rows_per_strip
        s0 = row0 // rps
        s1 = (row1 - 1) // rps
        rows_out = []
        with open(self.path, "rb") as f:
            for s in range(s0, s1 + 1):
                f.seek(int(self.strip_offsets[s]))
                buf = f.read(int(self.strip_byte_counts[s]))
                strip_rows = min(rps, self.length - s * rps)
                arr = np.frombuffer(buf, dtype=self.sample_dtype).reshape(
                    strip_rows, self.width, 2)
                rs = max(0, row0 - s * rps)
                re = min(strip_rows, row1 - s * rps)
                if re > rs:
                    rows_out.append(arr[rs:re, col0:col1, :])
        block = np.concatenate(rows_out, axis=0)  # (R, C, 2)：I,Q
        complex = (block[:, :, 0].astype(np.float32)
                   + 1j * block[:, :, 1].astype(np.float32))
        return complex.astype(np.complex64)


# --------------------------------------------------------------------------
# 4. 可选 ESA 辐射定标（ratio-to-reference）
# --------------------------------------------------------------------------
def apply_calibration(slc, cal_xml, cal_type, row0, col0):
    """按 ESA 比值法把 DN 复读数定标为 gamma0/sigma0/beta0 的复数：
        calibrated = DN * sqrt(LUT / |dn|^2)，使 |calibrated|^2 ≈ 所选 LUT。
    """
    import xml.etree.ElementTree as ET
    from scipy.interpolate import RegularGridInterpolator

    r = ET.parse(cal_xml).getroot()
    vecs = list(r.iter("calibrationVector"))
    lines, pixels, lut_grid, dn_mag2 = [], [], [], []
    for v in vecs:
        lines.append(int(v.findtext("line")))
        pix = np.array([int(x) for x in (v.findtext("pixel") or "").split()], dtype=np.int64)
        lut = np.array([float(x) for x in (v.findtext(cal_type) or "").split()], dtype=np.float64)
        dn_txt = (v.findtext("dn") or "0 0").split()
        dn = complex(float(dn_txt[0]), float(dn_txt[1]) if len(dn_txt) > 1 else 0.0)
        pixels.append(pix)
        lut_grid.append(lut)
        dn_mag2.append(abs(dn) ** 2)

    lines = np.array(lines, dtype=np.float64)
    pixels = pixels[0]  # 各 vector 的 pixel 网格一致
    lut_grid = np.array(lut_grid, dtype=np.float64)          # (M, N)
    dn_mag2 = np.array(dn_mag2, dtype=np.float64)            # (M,)

    interp = RegularGridInterpolator(
        (lines, pixels), lut_grid, bounds_error=False, fill_value=None)
    R, C = slc.shape
    RR, CC = np.meshgrid(row0 + np.arange(R), col0 + np.arange(C), indexing="ij")
    pts = np.stack([RR.ravel(), CC.ravel()], axis=1)
    lut_blk = interp(pts).reshape(R, C)

    dn2_blk = np.interp(row0 + np.arange(R).astype(float), lines, dn_mag2)
    dn2_blk = np.broadcast_to(dn2_blk[:, None], (R, C))

    factor = np.sqrt(lut_blk / dn2_blk).astype(np.complex64)
    return slc * factor


# --------------------------------------------------------------------------
# 5. ROI / 参数 / 输出
# --------------------------------------------------------------------------
def extract_rois(slc, roi_size, num_rois):
    R, C = slc.shape
    hh = min(int(roi_size[0]), R)
    ww = min(int(roi_size[1]), C)
    rois = []
    for k in range(num_rois):
        t = 0.5 if num_rois == 1 else (k + 0.5) / num_rois
        r0 = max(0, min(R - hh, int((R - hh) * t)))
        c0 = max(0, min(C - ww, int((C - ww) * t)))
        rois.append((r0, c0, hh, ww))
    return rois


def build_params(px, safe_dir, pol, shape, win, calibrate, cal_type):
    params = dict(px)
    params["product_id"] = os.path.basename(safe_dir.rstrip("/\\"))
    params["sensor"] = px.get("missionId", "Sentinel-1")
    params["mode"] = px.get("mode")
    params["swath"] = px.get("swath")
    params["product_type"] = px.get("productType")
    params["polarization"] = pol.upper()
    params["image_rows"] = int(shape[0])
    params["image_cols"] = int(shape[1])
    params["read_window"] = {"row0": int(win[0]), "row1": int(win[1]),
                             "col0": int(win[2]), "col1": int(win[3])}
    params["full_rows"] = px.get("numberOfLines")
    params["full_cols"] = px.get("numberOfSamples")
    params["calibrated"] = bool(calibrate)
    params["calibration_type"] = cal_type if calibrate else None
    params["generated_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    params["note"] = "第2周交付：SLC 复数读取 + 参数保存 + ROI；可选 ESA 辐射定标"
    return params


def _jsonable(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.complexfloating):
        return {"re": float(o.real), "im": float(o.imag)}
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(x) for x in o]
    return o


def write_json(obj, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(_jsonable(obj), f, ensure_ascii=False, indent=2)


def save_preview(slc, path, db_range=60):
    from PIL import Image

    amp = np.abs(slc)
    amp = np.where(amp > 0, 20 * np.log10(amp + 1e-12), -db_range)
    vmax = float(np.nanmax(amp))
    vmin = vmax - db_range
    img = np.clip((amp - vmin) / (vmax - vmin + 1e-12), 0, 1) * 255
    Image.fromarray(img.astype(np.uint8)).save(path)


# --------------------------------------------------------------------------
# 6. 主入口
# --------------------------------------------------------------------------
def read_slc(safe_dir, polarization="VV", out_dir=None,
             row0=0, row1=None, col0=0, col1=None,
             calibrate=False, cal_type="gamma",
             save=True, save_mat=False, plot=True,
             roi_size=(1024, 1024), num_rois=3,
             root_dir="D:/study/SAR"):
    ann, meas, cal = find_pol_files(safe_dir, polarization)
    params_xml = parse_annotation(ann)
    tiff = SlcTiff(meas)
    if row1 is None:
        row1 = tiff.length
    if col1 is None:
        col1 = tiff.width

    slc = tiff.read_block(row0, row1, col0, col1)

    if calibrate:
        if cal is None:
            raise RuntimeError("请求辐射定标但未找到 calibration XML")
        slc = apply_calibration(slc, cal, cal_type, row0, col0)
        params_xml["calibrated"] = True
        params_xml["calibration_type"] = cal_type

    params = build_params(params_xml, safe_dir, polarization, slc.shape,
                          (row0, row1, col0, col1), calibrate, cal_type)
    rois = extract_rois(slc, roi_size, num_rois)
    params["rois"] = [{"row": int(r[0]), "col": int(r[1]),
                       "rows": int(r[2]), "cols": int(r[3])} for r in rois]

    if save:
        out_dir = out_dir or os.path.join(root_dir, "data", "processed", "slc")
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(os.path.join(out_dir, "roi"), exist_ok=True)
        pol = polarization.upper()
        tag = f"{pol}_{slc.shape[0]}x{slc.shape[1]}"
        npy = os.path.join(out_dir, f"slc_{tag}.npy")
        np.save(npy, slc)
        jpath = os.path.join(out_dir, f"slc_{tag}_params.json")
        write_json(params, jpath)

        for k, (r0i, c0i, hh, ww) in enumerate(rois):
            sub = slc[r0i:r0i + hh, c0i:c0i + ww]
            np.save(os.path.join(out_dir, "roi", f"roi_{k}.npy"), sub)

        if save_mat:
            try:
                from scipy.io import savemat
                mpath = os.path.join(out_dir, f"slc_{tag}.mat")
                # 仅存复数图像；参数以 JSON 为准（避免 None 等不可序列化值）
                savemat(mpath, {"slc": slc})
            except Exception as e:  # pragma: no cover
                print("[warn] 保存 .mat 失败：", e)

        if plot:
            try:
                save_preview(slc, os.path.join(out_dir, f"slc_{tag}_preview.png"))
            except Exception as e:  # pragma: no cover
                print("[warn] 预览 PNG 失败：", e)

        params["output_npy"] = npy
        params["output_params_json"] = jpath
        write_json(params, jpath)  # 回填输出路径

    return slc, params


# --------------------------------------------------------------------------
# 7. CLI
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Sentinel-1 SLC 读取器（A 组 第2周）")
    ap.add_argument("--safe_dir", required=True, help="SLC 产品 .SAFE 目录")
    ap.add_argument("--polarization", default="VV")
    ap.add_argument("--row0", type=int, default=6000)
    ap.add_argument("--row1", type=int, default=None)
    ap.add_argument("--col0", type=int, default=6000)
    ap.add_argument("--col1", type=int, default=None)
    ap.add_argument("--calibrate", action="store_true", help="ESA 辐射定标")
    ap.add_argument("--cal_type", default="gamma", choices=["gamma", "sigma", "beta"])
    ap.add_argument("--out_dir", default=None)
    ap.add_argument("--save_mat", action="store_true")
    ap.add_argument("--no_plot", action="store_true")
    ap.add_argument("--roi_size", type=int, default=[512, 512], nargs=2, metavar=("H", "W"))
    ap.add_argument("--num_rois", type=int, default=3)
    ap.add_argument("--root_dir", default="D:/study/SAR")
    args = ap.parse_args()

    slc, params = read_slc(
        args.safe_dir, polarization=args.polarization,
        row0=args.row0, row1=args.row1, col0=args.col0, col1=args.col1,
        calibrate=args.calibrate, cal_type=args.cal_type,
        save=True, save_mat=args.save_mat, plot=not args.no_plot,
        roi_size=tuple(args.roi_size), num_rois=args.num_rois,
        root_dir=args.root_dir)

    print("SLC 读取完成：", slc.shape, "| 保存至", params.get("output_params_json"))
    print("fc=%.6g Hz, λ=%.5f m, PRF=%.3f Hz, 行=%d 列=%d" % (
        params.get("carrier_frequency_hz", 0), params.get("wavelength_m", 0),
        params.get("prf_hz", 0), params.get("image_rows", 0),
        params.get("image_cols", 0)))


if __name__ == "__main__":
    main()
