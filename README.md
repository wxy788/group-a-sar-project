# 哨兵 SAR 图像处理项目（A 组）

> 人工智能课程本科生大作业 · 围绕 Sentinel‑1 条带模式（Stripmap, SM）数据的 SAR 成像、深度自聚焦与深度旁瓣抑制。

## 1. 项目简介

本项目以 vibe coding 为主要开发方式，构建一套支持「数据导入 → 算法运行 → 结果对比 → 指标计算 → 导出」的 SAR 图像处理软件。仅处理 Sentinel‑1 SM 数据，最终形成可运行、可验收的 V1.0。

三个核心模块：
- **成像与数据**（MATLAB）：RD 成像、SM SLC 读取与参数保存。
- **深度自聚焦**（Python/PyTorch）：选定一种深度网络，恢复散焦复数图像，与 PGA 对比。
- **深度旁瓣抑制**（Python/PyTorch）：选定一种深度网络，抑制旁瓣并保留主瓣，与 SVA 切换。

## 2. 环境依赖

- MATLAB（含 **MATLAB Engine for Python**，本项目当前未装）
- Python 双 venv（均 uv 管理、CPU 版 torch，路径已被 .gitignore 排除）：
  - **`.venv`（Python 3.12）** → SVNN 旁瓣抑制 + 通用脚本（`read_slc.py` 等）
  - **`envs/pafnet`（Python 3.11）** → PAFnet 自聚焦（`torchbox 1.1.27` + `torchsar 1.1.4` 手动装）
- 界面：Streamlit（第8周，经 MATLAB Engine 调用传统函数）
- 详见 `requirements.txt`（含两条线一键复现命令、torchsar 手动安装说明）；
  自聚焦/旁瓣抑制各自的官方仓库已克隆到 `third_party/`。

## 3. 目录结构

```
SAR-Project/
├── README.md                 # 本文件
├── .gitignore
├── requirements.txt          # Python 依赖
├── docs/                     # 文档
│   ├── requirements.md       # 需求清单
│   ├── weekly-reports/       # 周报
│   └── references/           # 文献与笔记
├── data/
│   ├── raw/                  # SM SLC 原始产品（体积大，勿提交）
│   └── processed/            # 裁块、退化样本、参数 json
├── src/
│   ├── imaging/              # MATLAB: RD 成像、PGA、SVA、Level-0 解码
│   │   ├── read_level0.m     #   【本项目封装】Level-0 RAW 解码+复数读取+参数/ROI 导出
│   │   ├── read_slc.py       #   【本项目封装】Level-1 SLC 复数分块读取+定标+ROI
│   │   └── sentinel_l0_toolkit/  # 助教提供的解码/成像工具包（.venv 仍在 data/raw/Sentinel-L0条带）
│   ├── autofocus/            # Python/PyTorch: 自聚焦网络
│   ├── sidelobe/             # Python/PyTorch: 旁瓣抑制网络
│   └── gui/                  # Streamlit 界面 + MATLAB Engine 调用
├── envs/                     # Python 虚拟环境（.venv 已忽略；此目录见 .gitignore）
│   └── pafnet/               #   PAFnet 自聚焦环境（Python 3.11）
├── third_party/             # 官方仓库克隆（已忽略）：AutofocusSAR(PAFnet)、SVNN
├── models/                   # 训练好的权重与配置
├── results/
│   ├── figures/              # PNG 预览
│   └── metrics/              # CSV 指标
├── scripts/                  # 启动脚本、数据处理脚本
└── tests/                    # 复现 / 独立测试
```

## 4. 使用流程

1. 选择「仿真回波成像」或「SM SLC 导入」入口，查看产品参数。
2. 选择传统或网络方法，设置处理参数。
3. 运行算法，比较图像与剖面，计算指标并导出结果。
4. 各算法可独立运行，也可按「成像/SLC 导入 → 自聚焦 → 旁瓣抑制」串联。

## 5. 结果保存约定

- 复数结果：MAT；预览：PNG；指标：CSV；参数及处理链：JSON。
- 保存模型配置、权重、数据划分清单及必要日志（输入标识、算法版本、处理顺序、耗时）。
- 使用 Git 管理版本；训练与推理脚本随软件一并交付。

## 6. 子模块说明

| 模块 | 选型（已定稿） | 参考 |
|---|---|---|
| 自聚焦 | PAFnet（RAFnet 备选） | docs/references 中任务二笔记 |
| 旁瓣抑制 | SVNN（CNN+RSVA 备选） | docs/references 中任务三笔记 |
| 成像与数据 | sentinel_l0_toolkit（助教）+ read_level0.m（本项目封装） | src/imaging/README.md |

> **第 2 周数据线已完成**：环境部署（`requirements.txt`）、读取代码、以及
> 「按场景固定训练/验证/测试划分 + 标签与网络输入核对」详见
> [docs/data-pipeline.md](docs/data-pipeline.md)。
> 划分清单：`data/splits/scene_inventory.csv` → `data/splits/scene_split.json`；
> 核对脚本：`scripts/make_scene_split.py`、`scripts/verify_io_and_split.py`。

## 7. 验收要点

- 三模块均可运行；自聚焦与旁瓣抑制各完成一种深度网络（含训练/微调记录、可加载权重、独立测试、软件接口）。
- 真实 SM SLC 适配实验必须提交；无真值时分析图像变化与失败案例。

## 8. 参考资料

- 数据：[1] Sentinel‑1 Products；[2] Sentinel‑1 Processing；[3] Copernicus Browser。
- 自聚焦：[6] PAFnet (Liu et al. 2022, TGRS)；[7] RAFnet (Wu et al. 2026, Sensors)。
- 旁瓣抑制：[8] SVNN (Suo et al. 2024, TGRS)；[9] CNN+RSVA (Yuan et al. 2021, GRSL)。
