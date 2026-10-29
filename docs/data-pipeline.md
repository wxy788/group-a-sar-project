# 数据划分与 I/O 契约（第 2 周交付）

> 范围：第二周任务“按场景固定训练、验证和测试划分，核对标签与网络输入”。
> 配套脚本：`scripts/make_scene_split.py`、`scripts/verify_io_and_split.py`。
> 产物：`data/splits/scene_inventory.csv`（登记表）、`data/splits/scene_split.json`（权威划分）、`data/splits/verify_report_*.json`（核对报告）。

---

## 1. 按场景固定划分策略

### 为什么必须“按场景”而不是“按块/按像素”
两个网络都在**整景 SAR 图像**上提取训练样本（PAFnet 取 256×256 散焦块；SVNN 取沿方位/距离的 1D 剖面）。
若把这些块/剖面在 train/val/test 间**随机打散**，同一景的内容会同时出现在多个集合，
造成**数据泄漏**，测试指标虚高。因此划分的最小单元必须是 **scene（一次独立采集）**，
同一景产生的全部样本只能落在同一个集合。这正是 PAFnet 官方 `data.yaml` 的做法
（`trainid/validid/testid` 直接按场景文件编号分配）。

### 划分规则（确定性 + 可追加）
- **数据源**：`data/splits/scene_inventory.csv` —— 每个 scene 一行，含 `assigned_split` 列。
- **显式优先**：`assigned_split` 已填 `train/val/test` 的，生成器原样采用。
- **哈希补全**：留空的 scene 由 `sha256(seed:scene_id)` 按固定比例（默认 70/15/15）桶分配，
  结果完全可复现。
- **append-only**：新增场景时，旧 scene 的分配**永不改动**，只在尾部追加，从根上杜绝泄漏。
- **当前状态**：仅有 1 景真实数据 `S1A_20250303T140022_072E7E`（RAW+SLC 同一次采集），
  已显式保留为 **`test`**（两网络的真实泛化/定性评估集，不参与合成训练）。
  `train`/`val` 待后续经 `scripts/download_s1_data.py` 下载更多场景或合成场景后补充。

### 复现命令
```bash
# 由 inventory 生成 scene_split.json（纯标准库，任意 Python 3.8+）
python scripts/make_scene_split.py

# 仅做完整性/泄漏检查（不写文件）
python scripts/make_scene_split.py --check

# 加景后：在 scene_inventory.csv 追加一行，重新生成即可（旧分配不变）
```

---

## 2. 标签与网络输入契约

### 2.1 SVNN（旁瓣抑制，监督）—— 来自 `third_party/SVNN`
| 项 | 定义 |
|---|---|
| 训练输入 `sample` | 带旁瓣的 1D 幅度剖面（长度 512），`.mat` 键 `sample` |
| 监督标签 `gt` | 同长度、旁瓣被抑制的“干净”幅度剖面，`.mat` 键 `gt` |
| 数据布局 | `data/train/sample/*.mat` + `data/train/gt/*.mat`（val 同理） |
| 网络张量 | `(B, 1, L)`，1D 卷积；输入/输出**长度守恒** |
| 推理（真实 SLC） | 复数 `signal`(Na,Nr)，按 512 长度逐行（方位/距离分别）处理，
  实/虚部分别过网再合成复数输出（`inference.py`） |
| 损失 | `trainLoss`：误差项 `|pre-gt|` + 结构项 + 抑制项（抑制 `sam>gt` 区域） |
| 标签语义核对 | `sample` 旁瓣能量 **>** `gt` 旁瓣能量（输入更“脏”，标签更“净”） |

> 注：仓库 `dataset.py` 的 `BasicDataset` 直接 `net(sample)`，要求喂入 `(B,1,L)`；
> 本项目数据加载器需把每张 `.mat` 的 1D 剖面整理成 `(B,1,L)`，与仓库对接。

### 2.2 PAFnet（自聚焦，无监督）—— 来自 `third_party/AutofocusSAR/PAFnet`
| 项 | 定义 |
|---|---|
| 训练输入 | 散焦复数块 `(Na, Nr, 2)`（实部/虚部双通道），块尺寸 **256×256** |
| 数据文件 | `…_Defocused_Vr(…)_8000x256x256.h5`，键 `SI`(复数)、`ca`(相位系数)、`cr` |
| 训练监督 | **无**：以最小图像熵（Entropy）为损失，自聚焦无需成对标签 |
| 评估标签 `ca` | 真实方位相位误差多项式系数（如 `[2,3]` 次），仅供评估/可视化，
  **不参与训练损失** |
| 网络张量 | 输入 `(N,256,256,2)` → 输出 `(N,256,256,2)`，形状守恒；
  估计的 `ca` 形状 `(N, Ma)` |
| 真实数据说明 | 本项目的 Sentinel-1 真实 SLC 无 `ca` 真值，评估改用熵/对比度等指标 |

---

## 3. 核对结果（可运行验证）

`scripts/verify_io_and_split.py` 三模式，均写出 JSON 报告到 `data/splits/`：

```bash
python scripts/verify_io_and_split.py --net split    # 泄漏检查（纯标准库）
.venv/Scripts/python.exe        scripts/verify_io_and_split.py --net svnn   # SVNN 契约
envs/pafnet/Scripts/python.exe  scripts/verify_io_and_split.py --net pafnet # PAFnet 契约
```

最近一次核对（2026-10-29）：
- `split`：✅ 无跨集泄漏；⚠️ train/val 为空（待补充场景，符合设计）。
- `svnn`：✅ 5/5 —— 形状一致、float 幅度、前向保长、复数推理守恒、标签语义（sample 旁瓣 > gt）。
- `pafnet`：✅ 5/5 —— 输入 `(1,256,256,2)`、输出守恒、`ca(1,2)` 正确、熵为有限标量、聚焦后熵↓。

---

## 4. 后续衔接（第 3–4 周）
1. **造训练样本**（服务第 4 周）：
   - 自聚焦：从 SLC/RAW 抽 ROI，注入可控方位相位误差 → 散焦块 (`src/autofocus/`)。
   - 旁瓣抑制：抽强/弱目标 1D 剖面，合成 `sample`(高旁瓣)/`gt`(干净) 对 (`src/sidelobe/`)。
2. **样本 → 集合映射**：每个样本记录其 `scene_id`，由 `scene_split.json` 决定落入
   train/val/test，保证场景级不跨集（脚本已支持，待样本就绪后批量标注）。
3. **指标**：传统 PGA/SVA 与两网络的 PSLR/ISLR/熵对比（`results/metrics/`）。
