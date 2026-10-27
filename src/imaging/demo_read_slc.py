# -*- coding: utf-8 -*-
"""demo_read_slc.py — read_slc 用法示例（A 组 第2周 SLC 读取线）

运行方式：
  1) 命令行： python demo_read_slc.py
  2) 用解码 venv 的 python 直接跑（已含 numpy/scipy/tifffile/Pillow）：
     data/raw/Sentinel-L0条带/.venv/Scripts/python.exe demo_read_slc.py

说明：读取 VV 极化一个 2048×2048 块（避开图像零边界），保存
      复数 .npy + 参数 .json + 3 个 ROI .npy + 幅度预览 PNG（+可选 .mat）。
"""
import numpy as np
from read_slc import read_slc

SAFE = r"D:/study/SAR/data/raw/S1A_S2_SLC__1SDV_20250303T140022_20250303T140047_058141_072E7E_5C22.SAFE"

if __name__ == "__main__":
    slc, params = read_slc(
        SAFE, polarization="VV",
        row0=6000, row1=8048,          # 2048 行
        col0=6000, col1=8048,          # 2048 列
        calibrate=False,
        save=True, save_mat=True, plot=True,
        roi_size=(512, 512), num_rois=3,
    )
    print("done:", slc.shape, "complex64, 幅度均值=%.3f" % float(np.abs(slc).mean()))
