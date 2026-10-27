# src/imaging · 成像与数据模块

本目录存放**成像与数据**相关的代码（规划文档模块一：SAR 成像、SM SLC 读取、Level-0 解码、PGA/SVA 传统算法）。

## 目录结构

```
src/imaging/
├── read_level0.m                  # 【本项目封装】Level-0 RAW 解码 + 复数读取 + 参数/ROI 导出
├── sentinel_l0_toolkit/            # 助教提供的解码/成像工具包（原在 data/raw/Sentinel-L0条带）
│   ├── read_s1_raw.m               #   主入口：调 Python 解码器，输出复数回波 MAT + info
│   ├── IQ_correction.m             #   I/Q 校正（按 README 顺序在 SWST 对齐之前）
│   ├── echo_SWST.m                 #   按 SWST 对齐距离采样窗
│   ├── Read_Data_RD.m              #   RD 成像（距离压缩→RCMC→方位压缩），第3周用
│   ├── Shift_Delay_correction.m
│   ├── doppler_geo_parameter.m
│   ├── interpolate_s1_orbit.m
│   ├── decode_s1_raw.py            #   Python 解码器（被 read_s1_raw 经 system() 调用）
│   ├── decode_s1_navigation.py
│   ├── 设置SAR系统参数和成像参数.txt  # SARPara 结构体模板（供 Read_Data_RD 使用）
│   ├── requirements.txt            #   sentinel1decoder==2.1.0 等
│   └── README.md                   #   助教原版说明（权威，以它为准）
└── README.md                       # 本文件
```

> 工具包内 **`.venv` 仍保留在 `data/raw/Sentinel-L0条带/.venv`**（环境，已被 .gitignore 排除）。`read_level0.m` 会自动定位它；如移动过 venv，请同步修改 `locateVenv()` 的候选路径。

## 解码环境（Python venv）

`read_s1_raw.m` 经 `system()` 调用同目录 `decode_s1_raw.py`，需要 **Python 3.12** 环境
（依赖 `sentinel1decoder==2.1.0`、`numpy`、`scipy`、`h5py`、`hdf5storage`，见工具包 `requirements.txt`）。
该环境固定为 `data/raw/Sentinel-L0条带/.venv`（不进仓库）。

**本机已修复（2026-10-27）**：原 venv 的 `pyvenv.cfg` 指向不存在的 `codex-runtimes` base，
曾导致 `RAW decoding failed`。已把 `home`/`executable` 重指向本机 uv 管理的 CPython 3.12.14：

```
C:\Users\www\AppData\Roaming\uv\python\cpython-3.12-windows-x86_64-none\python.exe
# 前端：c:/Users/www/.local/bin/python3.12.exe（版本 3.12.14，与 venv 完全一致）
```

版本一致 → 已装的 cp312 编译包（scipy/h5py 等）ABI 兼容，直接可用，无需重装。

若换机或 venv 失效，重建方式：

```powershell
# 用真实 Python 3.12 作为 base 重建
<python3.12> -m venv data/raw/Sentinel-L0条带/.venv
data/raw/Sentinel-L0条带/.venv/Scripts/python.exe -m pip install -r src/imaging/sentinel_l0_toolkit/requirements.txt
```

> `read_level0.locateVenv()` 优先找 `sentinel_l0_toolkit/.venv`，其次 `data/raw/Sentinel-L0条带/.venv`；
> 移动 venv 后同步改候选路径。

## 数据位置（不进仓库，data/raw/** 已忽略）

| 数据 | 路径 |
|---|---|
| Level-0 RAW 产品 | `data/raw/S1A_S2_RAW__0SDV_20250303T140022_20250303T140047_058141_072E7E_7DD0/` 内含 `*SAFE/s1a-s2-raw-s-{vv,vh}-...dat` |
| Level-1 SLC 产品 | `data/raw/S1A_S2_SLC__1SDV_20250303T140022_20250303T140047_058141_072E7E_5C22.SAFE/` |
| Python 解码环境 | `data/raw/Sentinel-L0条带/.venv` |
| 处理输出 | `data/processed/level0/`（原始解码、校正块、参数 JSON、ROI） |

> 本产品为 **SDV = VV + VH 双极化**，与助教 README 示例（SDH=HH/HV）不同——调用时务必指定 `Polarization`（见下）。

## 运行 Level-0 读取（read_level0.m）

在 **MATLAB** 中，将当前目录切到 `src/imaging`（或确保它在路径上），然后：

```matlab
% 1) 读一段 4000 行块（含 IQ 校正 + SWST 对齐 + 3 个 ROI + 参数 JSON）
read_level0('RawRoot', ...
  'D:/study/SAR/data/raw/S1A_S2_RAW__0SDV_20250303T140022_20250303T140047_058141_072E7E_7DD0', ...
  'Polarization', 'VV');

% 2) 整景磁盘模式解码（~10GB，需 LoadIntoMemory=false；不做 IQ 校正）
read_level0('RawRoot', 'D:/study/SAR/data/raw/S1A_S2_RAW__0SDV_20250303T140022_20250303T140047_058141_072E7E_7DD0', ...
  'Polarization', 'VV', 'Full', true);

% 3) 自定义 ROI（行/列/高/宽，相对读出的块坐标）
read_level0('RawRoot', '...', 'Polarization','VV', 'ROIs', {[200 100 1024 1024], [2000 200 1024 1024], [4000 300 1024 1024]});
```

输出：
- `data/processed/level0/level0_vv_start1_n4000_raw.mat`：read_s1_raw 原始解码（含 `iq`、`timeline`、`radar`、`orbit` 等）
- `..._corrected.mat`：IQ 校正 + SWST 对齐后的块
- `..._params.json`：产品标识、模式、极化、fc/lambda/PRI/PRF/采样率、`Na_all`、ROI 位置等参数
- `roi/roi_k.mat` + `roi_k.png`：3 个 ROI

## 关键约定（务必遵守）

- **极化**：本产品 SDV → 传 `'VV'`（或 `'VH'`）；不传默认 `'HH'` 会找不到文件。
- **内存**：单极化整景 `56435×24424` complex single ≈ 10.27 GiB；整景务必 `Full=true`（磁盘模式）。小块的 `NumPulses`（默认 4000）约 0.8 GiB，可载入内存做校正。
- **零行掩码**：标定/噪声/缺包 PRI 在 `iq` 中为全零行，用 `info.timeline.echo_mask` 区分，训练/成像时剔除。
- **不写回 SAFE**：输出一律在 `data/processed/`，绝不写入原始产品目录。
- **每调用一次都重新解码**：完整整景解码较慢，按需使用；已有 MAT 可经 `matfile()` 直接读块。

## SLC 读取（read_slc.py，第2周另一条线 · 已完成）

`read_slc.py` 读取 SLC 测量 TIFF 的**复数图像**（实测格式：I/Q 交织 int16 complex、4 字节/像素、
小端、无压缩、逐行 strip），解析 `annotation/*.xml` 成像参数，可选 ESA 辐射定标，
抽取 ROI 训练样本，输出复数 `.npy` + 参数 `.json`（+可选 MATLAB `.mat`、幅度预览 PNG）。
按行 strip **按需读取**，整景（39297×20104 ≈ 6.3 GB complex64）不必载入内存。

```bash
# 命令行（用解码 venv 的 python，已含 numpy/scipy/tifffile/Pillow）
data/raw/Sentinel-L0条带/.venv/Scripts/python.exe src/imaging/read_slc.py \
  --safe_dir "data/raw/S1A_S2_SLC__1SDV_20250303T140022_20250303T140047_058141_072E7E_5C22.SAFE" \
  --polarization VV --row0 6000 --row1 8048 --col0 6000 --col1 8048
```

```python
from read_slc import read_slc
slc, params = read_slc(SAFE, 'VV', row0=6000, row1=8048, col0=6000, col1=8048)  # complex64
```

输出（`data/processed/slc/`）：`slc_<POL>_<H>x<W>.npy / .mat`、`..._params.json`、
`..._preview.png`、`roi/roi_k.npy`（训练/测试样本块）。

要点：
- **极化**：SDV 产品，VV/VH 各一套 annotation + measurement + calibration XML。
- **零边界**：图像边缘存在零行/零列，取块时避开（demo 从 6000 起）。
- **定标**：`calibrate=True` 按 ESA 比值法 `DN·sqrt(LUT/|dn|²)` 得复数，其模方 ≈ gamma0/sigma0/beta0
  （默认关闭；自聚焦/旁瓣抑制任务用原始复数 DN 即可）。
- **参数 JSON**：fc/λ/PRF/像素间距/采样率/入射角/平台航向/UTC 时间 + **annotation 内 17 个轨道状态向量**。
- 已实测：2048×2048 VV 块（6000:8048），幅度均值 ≈ 40.4，预览为清晰地形 SAR 图像；定标路径已验证。
- 依赖：numpy（必需）；tifffile、scipy、Pillow（已在解码 venv 中，见根 `requirements.txt`）。
