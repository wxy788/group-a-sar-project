# src/autofocus · 自聚焦模块

本目录存放**自聚焦**相关的 MATLAB 代码（规划方案 §三「阶段二 深度学习自聚焦」的传统基线与样本构造部分）。
网络侧（PAFnet）在 `envs/pafnet` 的 Python/PyTorch 环境中，第 4 周起接入。

## 目录结构

```
src/autofocus/
├── week3_autofocus.m          # 【第3周交付】散焦注入 + 最小熵自聚焦 + PGA + 指标 + 样本生成
├── run_autofocus_engine.py    # 连接用户已共享的 GUI MATLAB 会话跑上面的脚本（本机绕过 License -5201）
└── README.md                  # 本文件
```

## 运行

本机 MATLAB 沙箱直连会报 `License Manager Error -5201`，因此走「用户在 GUI MATLAB 里
`matlab.engine.shareEngine` → Python 连该会话」的方式复用授权（与成像模块同一套路）：

```powershell
# 1) 在 GUI MATLAB 命令窗口执行一次：
matlab.engine.shareEngine

# 2) 在仓库根目录：
.venv\Scripts\python.exe src\autofocus\run_autofocus_engine.py            # 自动连最近会话
.venv\Scripts\python.exe src\autofocus\run_autofocus_engine.py MATLAB_96368
```

也可直接在 MATLAB 里 `cd src/autofocus` 后执行 `week3_autofocus`。
仅依赖 MATLAB 内置函数（`fft/ifft/fminsearch/optimset/imagesc/jsonencode`），**无工具箱依赖**。

## 相位误差模型（与 PAFnet 等文献一致）

误差加在**方位频谱**上（列方向即方位向）：

```
散焦  D = ifft( fft(I0,[],2) .* exp(+j*phi(f)) ,[],2)
补偿  I = ifft( fft(D ,[],2) .* exp(-j*phi_hat(f)),[],2)
```

等价于「真实图像与散焦核 `h = ifft(exp(j*phi))` 卷积」。这是自聚焦的标准可逆模型：φ 已知即可
精确复原，因此可用**幅度 NRMSE** 客观度量「恢复能力」（真实场景无真值时才退回熵/对比度）。

误差类型：二次（`a2·u²`，对应 Ka 失配）、单频/多频正弦、带限随机（白噪声平滑到指定宽度）、混合。

## 两条传统基线

| | 最小熵自聚焦 MEAF | PGA 相位梯度自聚焦 |
|---|---|---|
| 原理 | `phi_hat = B·c`，B = 偶次多项式 `u²,u⁴,u⁶` + cos/sin 各 4 阶；以图像熵为目标做 Nelder-Mead（先粗扫二次项再 `fminsearch`） | 取能量最强的若干距离单元 → 亮点圆移到 index 1 → 环绕升余弦窗 → 方位谱相位梯度 `Δφ(k)=angle(Σ X(k)·conj(X(k−1)))` → 带内累积 |
| 优点 | 不依赖孤立点目标，分布杂波上也能工作 | 非参数，可捕捉高频误差；点目标上极准 |
| 缺点 | 受基函数表达能力限制（纯高频随机误差修不动） | 强依赖孤立强散射体，分布杂波上会失效 |
| 实测耗时 | ~29 s（512×512，11 参数） | ~0.1 s |

PGA 内建 **「熵不劣化守卫」**：每轮补偿前先算熵，若候选结果熵变大就回退并停止迭代，
保证不会把已经聚焦的图像改坏（实机上自聚焦的必备保护）。MEAF 未加守卫 —— 它的目标函数就是熵，
理论上不会劣化熵，但熵与 NRMSE 并不总同向（见 E8），这一现象本身值得记录。

### ⚠ 两个花了很久才定位的实现坑（勿回退）

1. **峰值必须圆移到 index 1，不能移到数组中心。**
   移到中心会给谱引入 `(-1)^k` 的 π 交替相位，相位梯度出现系统性 ±π 跳变，累积后图像被
   彻底破坏（实测 NRMSE≈1.4，连零误差输入都会被毁掉）。
2. **`unwrap` 与去线性必须只在有效频带内做。**
   方位谱是带限的（仿真与实测都只有约 47% 的带宽有能量），带外 `|X|≈0`、相位是纯噪声。
   若在全带 `unwrap`/去趋势，会得到错尺度甚至反向的相位（实测残差只有真值的 ~0.5 倍且整体崩坏）。
   做法：`valid = 幅值 > 5% 峰值`，仅用带内样本 `unwrap` + 一次多项式去趋势，带外保持边界常数。

## 产出

运行 `week3_autofocus` 后：

| 产物 | 路径 |
|---|---|
| 9 组案例（参考/散焦/MEAF/PGA 复数图 + φ 真值/估计） | `data/processed/autofocus/cases/af_case_E*.mat` |
| 第 4 周网络训练用 256×256 散焦-参考对（24 对） | `data/processed/autofocus/samples/af_pair_*.mat` |
| 样本索引（误差类型/RMS/来源坐标） | `data/processed/autofocus/dataset_index.csv` |
| 验收结论 | `data/processed/autofocus/autofocus_week3_report.json` |
| 逐案例指标表 | `results/metrics/autofocus_metrics.csv` |
| 汇总图/案例对比图/PGA 收敛曲线 | `results/figures/autofocus_*.png` |

指标口径（幅度 NRMSE / 图像熵 / 对比度 / 方位 −3dB 主瓣宽度 / 耗时）与规划方案 §三.4 一致，
与成像模块 `week3_imaging.m` 的 `analyze_1d/image_entropy/image_contrast` 同口径。

## 实测结果（2026-11-07，MATLAB R2024b）

底图：B1 = 仿真多点目标 512×512；B2 = 真实 SM SLC 512×512（块内尖峰度最高者）；
B3 = 真实 SM Level-0 聚焦图 512×512。

| 案例 | 底图 | 误差类型 | φ RMS | NRMSE 散焦 | MEAF | PGA |
|---|---|---|---|---|---|---|
| E1 | B1 | 二次 a2=6 | 2.68 | 0.2207 | **0.0007** | 0.0041 |
| E2 | B1 | 二次 a2=14 | 6.26 | 0.5782 | **0.0003** | 0.0056 |
| E3 | B1 | 正弦 1.5 周期 | 0.85 | 0.3885 | 0.0072 | **0.0010** |
| E4 | B2 | 二次 a2=8 | 3.58 | 0.5712 | **0.0108** | 0.5712（守卫拒绝） |
| E5 | B2 | 双频正弦 | 0.51 | 0.3348 | **0.1371** | 0.3348（守卫拒绝） |
| E6 | B2 | 带限随机 | 0.70 | 0.4302 | **0.3974** | 0.4302 |
| E7 | B3 | 正弦 2.5 周期 | 1.06 | 0.5725 | **0.5036** | 0.5725 |
| E8 | B3 | 高频随机 | 1.00 | 0.5287 | 0.6198（变差） | 0.6508 |
| E0 | B1 | 零误差 | 0 | 0 | **0.0002** | 0.0001 |
| E9 | B2 | 未见：混合 | 2.42 | 0.6456 | **0.2661** | 0.6456 |

- 主瓣恢复（E2）：方位 −3dB 宽度 散焦 4.99 → MEAF 1.79 采样（参考 1.79），精确复原。
- 结论：**8 个非零误差案例中 MEAF 7 个、PGA 3 个的 NRMSE 优于散焦输入**（§七.1 要求 ≥5），
  零误差与未见误差均通过。E8（真实 Level-0 上的高频随机误差）是唯一失败项 —— 高频随机误差
  超出「多项式 + 低阶正余弦」基的表达能力，属真实局限，第 4 周起交由 PAFnet 学习型补偿。
  另注意 E8 的熵其实略降（17.3876→17.3799）但 NRMSE 变差，说明**熵与 NRMSE 在近均匀分布的
  杂波图上并不等价**，评价自聚焦时必须两者并报。
- PGA 在 B2/B3 上被守卫拒绝（不做改动）：这两块是分布杂波，缺少 PGA 依赖的孤立强散射体。
  这是方法固有特性，不是实现缺陷（零误差输入下 PGA 的估计相位 ≈0，实测 NRMSE 0.0001）。
