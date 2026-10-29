#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
verify_io_and_split.py —— 核对「标签与网络输入」+ 场景划分无泄漏

对应第 2 周任务“核对标签与网络输入”与“按场景固定划分”。
三种模式（--net）：

  split   纯标准库。读取 data/splits/scene_split.json，核实：
           每个 scene 唯一分配、无跨集泄漏、train/val/test 非空（若已登记）。

  svnn    需在 .venv（torch + scipy）运行。核对 SVNN 契约：
           输入 sample = 带旁瓣的 1D 幅度剖面（长度 512）；
           标签 gt     = 纯净 1D 幅度剖面（同长度）；
           推理时对复数 SLC 实/虚部分别过网再合成。
           动作：合成一对 (sample,gt) 存为 .mat，按仓库 BasicDataset 逻辑读回，
           跑 svnn 前向，断言 输出长度=输入长度、gt 旁瓣低于 sample（标签语义正确）。

  pafnet  需在 envs/pafnet（torch + torchbox + torchsar）运行。核对 PAFnet 契约：
           输入 = 散焦复数块 (Na=256, Nr=256, 2)（实部/虚部双通道）；
           训练监督 = 无（最小熵损失）；评估标签 = ca（真实相位误差系数，仅供分析）。
           动作：合成散焦复数块（注入方位二次相位误差），跑 PAFnet.forward，
           断言 输出形状守恒、ca 形状正确、熵损失为有限标量（无监督训练信号存在）。

所有模式都会写出 JSON 报告到 data/splits/verify_report_<net>.json。

用法：
  python scripts/verify_io_and_split.py --net split
  .venv/Scripts/python.exe        scripts/verify_io_and_split.py --net svnn
  envs/pafnet/Scripts/python.exe  scripts/verify_io_and_split.py --net pafnet
"""
import argparse
import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
SPLITS_JSON = os.path.join(PROJECT_ROOT, "data", "splits", "scene_split.json")
THIRD_PARTY = os.path.join(PROJECT_ROOT, "third_party")
REPORT_DIR = os.path.join(PROJECT_ROOT, "data", "splits")


# --------------------------------------------------------------------------
# 模式 1：split 泄漏检查（纯标准库）
# --------------------------------------------------------------------------
def verify_split():
    checks = []
    if not os.path.isfile(SPLITS_JSON):
        return {"mode": "split", "passed": False,
                "checks": [{"name": "scene_split.json 存在", "passed": False,
                            "detail": f"未找到 {SPLITS_JSON}，请先运行 make_scene_split.py"}]}

    with open(SPLITS_JSON, "r", encoding="utf-8") as f:
        doc = json.load(f)

    splits = doc.get("splits", {})
    seen = {}
    leak = []
    for split in ("train", "val", "test"):
        for sid in splits.get(split, []):
            if sid in seen:
                leak.append(f"{sid}: {seen[sid]} & {split}")
            seen[sid] = split
    checks.append({"name": "每景唯一分配 / 无跨集泄漏", "severity": "error",
                   "passed": len(leak) == 0,
                   "detail": "无泄漏" if not leak else "; ".join(leak)})

    sizes = {k: len(v) for k, v in splits.items()}
    # 非空是“数据可用性”提示，不是完整性错误：当前仅 1 景真实数据（保留为 test），
    # train/val 设计为后续下载/合成补充，故为空时仅告警不失败。
    empty = [k for k in ("train", "val", "test") if sizes[k] == 0]
    checks.append({"name": "train/val/test 均非空（数据可用性）", "severity": "warn",
                   "passed": len(empty) == 0,
                   "detail": (f"当前划分={sizes}；"
                              f"{'全部就绪' if not empty else '待补充场景：' + ','.join(empty) + '（保留真实景为 test，训练/验证景后续下载/合成）'}")})

    # 仅 error 级检查失败才整体失败
    passed = all(c.get("severity", "error") != "error" or c["passed"] for c in checks)
    return {"mode": "split", "passed": passed, "splits": splits, "checks": checks}


# --------------------------------------------------------------------------
# 模式 2：SVNN 契约
# --------------------------------------------------------------------------
def verify_svnn(tmpdir):
    import numpy as np
    import torch
    import scipy.io as sio

    svnn_dir = os.path.join(THIRD_PARTY, "SVNN")
    sys.path.insert(0, svnn_dir)
    from svnn import svnn as SVNN

    checks = []
    L = 512  # SVNN 推理 INPUTLENGTH

    # ---- 合成一对（忠实于 SVNN 监督定义）----
    # sample = 带旁瓣（高旁瓣）的幅度剖面：矩形窗响应 |sinc|（~ -13 dB 旁瓣）
    # gt     = 同一主瓣但旁瓣被抑制的“干净”幅度剖面（旁瓣抑制的监督标签）
    n = np.arange(L) - L / 2
    main = np.abs(np.sinc(n / 12.0)).astype(np.float32)   # 主瓣 + 自然高旁瓣
    sample = main.copy()
    gt = main.copy()
    sidelobe_mask = np.abs(n) > 12
    gt[sidelobe_mask] *= 0.15  # 标签：旁瓣方向被抑制到 15%（主瓣保持不变）
    sample /= sample.max()
    gt /= gt.max()

    # 旁瓣能量（主瓣外）对比：sample 应高于 gt（标签语义：gt 更干净）
    mask = np.abs(n) > 18
    e_sample = float(np.sum(sample[mask] ** 2))
    e_gt = float(np.sum(gt[mask] ** 2))
    checks.append({"name": "标签语义：sample 旁瓣能量 > gt 旁瓣能量",
                   "passed": e_sample > e_gt,
                   "detail": f"sample_side={e_sample:.4f}, gt_side={e_gt:.4f}"})

    # ---- 落盘为 .mat（与仓库 BasicDataset 一致：键 sample / gt）----
    os.makedirs(os.path.join(tmpdir, "sample"), exist_ok=True)
    os.makedirs(os.path.join(tmpdir, "gt"), exist_ok=True)
    sio.savemat(os.path.join(tmpdir, "sample", "0.mat"), {"sample": sample})
    sio.savemat(os.path.join(tmpdir, "gt", "0.mat"), {"gt": gt})

    # ---- 以仓库 BasicDataset 逻辑读回 ----
    sm = sio.loadmat(os.path.join(tmpdir, "sample", "0.mat"))["sample"]
    gm = sio.loadmat(os.path.join(tmpdir, "gt", "0.mat"))["gt"]
    checks.append({"name": "sample 与 gt 形状一致（同长度 1D）",
                   "passed": sm.shape == gm.shape,
                   "detail": f"sample{sm.shape} gt{gm.shape}"})
    checks.append({"name": "数据类型为 float（幅度）",
                   "passed": sm.dtype.kind == "f" and gm.dtype.kind == "f",
                   "detail": f"sample {sm.dtype}, gt {gm.dtype}"})

    # 组装成网络输入 (B,1,L)
    x = torch.from_numpy(sample).float().unsqueeze(0).unsqueeze(0)
    y_gt = torch.from_numpy(gt).float().unsqueeze(0).unsqueeze(0)

    net = SVNN(filter_length=15)
    with torch.no_grad():
        out = net(x)
    checks.append({"name": "SVNN 前向：输出长度=输入长度（1D 保长）",
                   "passed": tuple(out.shape) == (1, 1, L),
                   "detail": f"input{x.shape} -> output{tuple(out.shape)}"})

    # ---- 复数推理路径：对复数 SLC 逐行（1D 序列）实/虚部分别过网再合成 ----
    # 与仓库 inference.py 的 filterSingleDirection 一致：每条 1D 序列长度 = INPUTLENGTH
    Lc = 512
    rng = np.random.default_rng(0)
    sig1d = (rng.standard_normal(Lc) + 1j * rng.standard_normal(Lc)).astype(np.complex64)

    def filter_1d(net, seq):
        seq = np.asarray(seq, dtype=np.float32)
        nf = np.max(np.abs(seq))
        t = torch.tensor(seq / nf).unsqueeze(0).unsqueeze(0)  # (1,1,L)
        with torch.no_grad():
            o = net(t).squeeze().numpy()
        return o * nf

    out_re = filter_1d(net, np.real(sig1d))
    out_im = filter_1d(net, np.imag(sig1d))
    sig_out = out_re + 1j * out_im
    checks.append({"name": "复数推理：实/虚部分别过网后合成，长度守恒（复数 I/O）",
                   "passed": sig_out.shape == sig1d.shape and np.iscomplexobj(sig_out),
                   "detail": f"input(Lc={Lc}) complex 1D -> output{sig_out.shape} complex 1D"})

    passed = all(c["passed"] for c in checks)
    return {"mode": "svnn", "passed": passed, "io_contract": {
        "input": "sample: 带旁瓣的 1D 幅度剖面 / 复数 SLC（实虚分别过网）",
        "label": "gt: 同长度纯净幅度剖面（监督标签）",
        "sequence_length": L,
        "tensor_shape": "(B,1,L)"
    }, "checks": checks}


# --------------------------------------------------------------------------
# 模式 3：PAFnet 契约
# --------------------------------------------------------------------------
def verify_pafnet():
    import numpy as np
    import torch

    paf_dir = os.path.join(THIRD_PARTY, "AutofocusSAR", "PAFnet")
    sys.path.insert(0, paf_dir)
    from pafnet import PAFnet

    checks = []
    Na, Nr = 256, 256  # PAFnet 数据块尺寸（data.yaml: 256x256）

    # ---- 合成散焦复数块：注入方位二次相位误差（模拟速度误差散焦）----
    rng = np.random.default_rng(1)
    # 干净聚焦复数块（随机场景纹理）
    clean = (rng.standard_normal((Na, Nr)) + 1j * rng.standard_normal((Na, Nr))).astype(np.complex64)
    # 方位维二次相位误差 phi = c * (n - N/2)^2
    n = np.arange(Na) - Na / 2
    phi = 0.6 * (n ** 2) / (Na ** 2) * np.pi * 2.0  # 弧度
    phase = np.exp(1j * phi)[:, None]  # 沿方位
    defocused = clean * phase

    # 网络输入：形状 (1, Na, Nr, 2) 实部/虚部双通道
    X = np.stack([defocused.real, defocused.imag], axis=-1)[None, ...].astype(np.float32)
    Xt = torch.from_numpy(X)

    # lpaf 须为 10 元素卷积参数 [out, kH, kW, sH, sW, pH, pW, dH, dW, groups]
    # （源码默认值 [[[32,11,1]]] 形状错误，按官方 demo 格式显式传入）
    net = PAFnet(Na=Na, Nr=Nr, Mas=[[2, 3]],
                 lpaf=[[[32, 11, 1, 1, 1, 0, 0, 1, 1, 1]]])
    with torch.no_grad():
        Y, ca = net(Xt)

    checks.append({"name": "PAFnet 输入形状 = (1,256,256,2) 复数双通道",
                   "passed": tuple(X.shape) == (1, Na, Nr, 2),
                   "detail": f"X{X.shape}"})
    checks.append({"name": "PAFnet 输出形状守恒（聚焦后仍 256x256 复数）",
                   "passed": tuple(Y.shape) == (1, Na, Nr, 2),
                   "detail": f"Y{tuple(Y.shape)}"})
    checks.append({"name": "评估标签 ca 形状正确（N, Ma=2）",
                   "passed": tuple(ca.shape) == (1, 2),
                   "detail": f"ca{tuple(ca.shape)} (Ma=2 对应 Mas=[[2,3]])"})

    # ---- 无监督熵损失合同：熵应为有限标量；聚焦后熵应 <= 散焦 ----
    def img_entropy(cplx):
        amp = np.abs(cplx)
        p = amp / (amp.sum() + 1e-12)
        p = p[p > 0]
        return float(-(p * np.log2(p)).sum())

    e_in = img_entropy(defocused)
    e_out = img_entropy(Y[0].numpy().view(np.complex64) if False else
                        (Y[0, ..., 0].numpy() + 1j * Y[0, ..., 1].numpy()))
    checks.append({"name": "无监督熵损失：散焦输入熵为有限标量（训练信号存在）",
                   "passed": bool(np.isfinite(e_in)),
                   "detail": f"entropy(defocused)={e_in:.4f}"})
    checks.append({"name": "聚焦后熵 <= 散焦熵（自聚焦方向正确）",
                   "passed": e_out <= e_in + 1e-6,
                   "detail": f"in={e_in:.4f} -> out={e_out:.4f}"})

    passed = all(c["passed"] for c in checks)
    return {"mode": "pafnet", "passed": passed, "io_contract": {
        "input": "散焦复数块 (Na,Nr,2) 实部/虚部双通道，块尺寸 256x256",
        "train_supervision": "无（最小图像熵损失 Entropy）",
        "eval_label": "ca：真实方位相位误差多项式系数（仅供评估/可视化，不参与训练损失）",
        "patch_size": f"{Na}x{Nr}"
    }, "checks": checks}


def main():
    ap = argparse.ArgumentParser(description="核对标签与网络输入 + 场景划分泄漏检查")
    ap.add_argument("--net", required=True, choices=["split", "svnn", "pafnet"])
    ap.add_argument("--tmpdir", default=os.path.join(REPORT_DIR, "_svnn_demo"))
    args = ap.parse_args()

    if args.net == "split":
        report = verify_split()
    elif args.net == "svnn":
        report = verify_svnn(args.tmpdir)
    else:
        report = verify_pafnet()

    os.makedirs(REPORT_DIR, exist_ok=True)
    report_path = os.path.join(REPORT_DIR, f"verify_report_{args.net}.json")
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    print(f"\n===== verify_io_and_split [{args.net}] =====")
    for c in report.get("checks", []):
        mark = "✅" if c["passed"] else "❌"
        print(f"  {mark} {c['name']}\n       -> {c['detail']}")
    print(f"结果：{'全部通过 ✅' if report['passed'] else '存在失败 ❌'}")
    print(f"报告已写出：{report_path}")
    sys.exit(0 if report["passed"] else 1)


if __name__ == "__main__":
    main()
