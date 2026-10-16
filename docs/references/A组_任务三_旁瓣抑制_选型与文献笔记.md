# A组 · 任务三（旁瓣抑制方向）· 选型与文献笔记 + 输入输出说明

---

## 一、任务背景

SAR 成像响应为 sinc 形，强散射目标旁瓣会遮蔽邻近弱目标、破坏结构特征。目标：抑制旁瓣同时**保留主瓣能量、分辨率与弱目标**，形成可与传统 **SVA（空间变迹）** 切换的处理模块。

---

## 二、候选方案对比

### 2.1 SVNN（推荐首选）

- **文献**：[8] Yuxi Suo, Kun Fu, Youming Wu, Qingbiao Meng, Tian Miao, Xian Sun. *Spatially Variant Filtering Network Based on Generalized Optimal Constraints for Sidelobe Suppression in SAR Images.* IEEE TGRS, 2024, vol. 62, Art no. 5228314. DOI: 10.1109/TGRS.2024.3498592
- **动机**：传统空间变滤波在「滤波设计约束 + 最小幅度约束 MAC」下求解解析滤波器，约束过窄导致抑制有限，且 MAC 无差别最小化幅度会**不可避免造成主瓣损失**。SVNN 在更宽松约束下学习更优的空间变滤波器，并保持主瓣。
- **两个核心模块**：
  - **SVFG（Spatially Variant Filter Generation）**：自适应生成随位置变化的滤波器（一维 U‑Net 式滤波器生成，需复现其序列长度与下采样约束）；
  - **PSC（Paralleled Shifted Convolution）**：转换信号格式，实现快速、并行的空间变滤波。
- **监督**：以**最小误差约束 MEC（minimum error constraint）**监督，主瓣保持更优（优于无监督 MAC）。
- **输入 / 输出**：
  - 输入：复数图像（实/虚或幅度相位，按作者材料核对标签含义）；
  - 输出：旁瓣抑制后**复数图像**（保留复数输入输出，检查相位变化）。
- **代码与数据**：官方仓库 `github.com/suoyuxi/SVNN`（含 `svnn.py`、`dataset.py`、`train.py`、`inference.py`、`eval.py`）。论文在 GaoFen‑3、TerraSAR‑X 卫星数据上验证。
- **优点**：官方代码与材料齐备、监督训练、SOTA 主瓣保持，复现门槛最低。

### 2.2 CNN + RSVA（备选）

- **文献**：[9] Sen Yuan, Ze Yu, Chunsheng Li, Shusen Wang. *A Novel SAR Sidelobe Suppression Method Based on CNN.* IEEE GRSL, 2021, 18(1):132–136. DOI: 10.1109/lgrs.2020.2968336
- **核心思想**：传统 **RSVA（鲁棒空间变迹）** 最有效，但相位因素会导致杂波目标**能量损失**；用 **CNN 补偿 RSVA 的主瓣能量损失**。
- **输入 / 输出**：
  - 输入：RSVA 输出（或复数图像）；
  - 输出：能量补偿后的图像（SNR 保持更优）。
- **代码**：作者机构资料，**需自行联系或实现**（无统一公开仓库）。
- **注意**：选择此法需自行搭建 RSVA 与 CNN 流程，工作量高于 SVNN。

### 2.3 对比表

| 维度 | SVNN（推荐） | CNN + RSVA |
|---|---|---|
| 年份 / 出处 | 2024 IEEE TGRS | 2021 IEEE GRSL |
| 方法类型 | 空间变滤波网络（生成式滤波器） | 传统 RSVA + CNN 能量补偿 |
| 主瓣保持 | MEC 监督，保持好 | CNN 补偿 RSVA 能量损失 |
| 官方代码 | ✅ 有 | ⚠️ 需自行实现 |
| 上手难度 | 低 | 中（需自己写 RSVA） |

---

## 三、本组选型建议

**首选 SVNN**：官方代码与训练/推理材料齐全、监督训练、主瓣保持 SOTA，第 4 周跑通把握大。
**CNN+RSVA 作为备选**：若组内希望对比「传统 SVA 家族 + 网络补偿」的路线，可后期采用；注意需自行实现 RSVA 与 CNN 串联。

---

## 四、输入输出说明（以 SVNN 为准，待定稿后替换）

- **输入**：`H×W` 复数图像 `I_in`（建议训练块 128×128 或 256×256，遵循网络序列长度约束）。
- **处理流程**：SVFG 生成随位置变化的滤波器 → PSC 并行空间变卷积 → 输出旁瓣抑制复数图像。
- **输出**：`H×W` 复数图像 `I_out`，保留复数；与传统 **SVA** 可切换对比。
- **物理机理**：旁瓣降低 vs 主瓣保真需平衡；重点检查弱目标是否保留、是否产生虚假峰、相位是否异常变化。

---
