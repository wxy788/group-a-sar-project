# 在 MATLAB 中读取 Sentinel-1 完整 RAW 回波

`read_s1_raw.m` 调用本目录的 Python 解码器，提取所选时间范围内的所有 ECHO acquisition chunks，按照包头 PRI Count 放回完整慢时间轴。成像区间定义为首个 ECHO 包到最后一个 ECHO 包，区间内的标定、噪声和缺包时隙均保留为零行，使用独立掩码区分。

## 使用

在 MATLAB 中将当前目录切换到本工具目录。默认读取本地洛杉矶 S1A S2 HH 文件的完整成像 PRI 时间轴：

```matlab
[echo, info] = read_s1_raw();
size(echo)                     % 56435 x 24424
info.image.Na_all              % 56435，整个成像区间的 PRI 数
% echo 是 single 复数矩阵，行=PRI 时隙，列=距离采样点。
```

单极化完整成像矩阵本身需要约 **10.27 GiB 内存**，加载及后续计算还需要额外空间。本机为 16 GB 内存，建议使用磁盘模式：它仍解码并保存全部回波，只按需加载矩阵片段。

```matlab
[data, info] = read_s1_raw('', Inf, 'LoadIntoMemory', false, 'Plot', false);
size(data, 'iq')               % 完整磁盘矩阵的尺寸
block = data.iq(1:2200, :);     % 按需读入一段，保留该段内零行

% 另一极化，同样提取所有回波块
[dataHV, infoHV] = read_s1_raw('', Inf, 'Polarization', 'HV', ...
    'LoadIntoMemory', false, 'Plot', false);

% 读取另一份 RAW 测量文件（不是 annot.dat / index.dat）
[data, info] = read_s1_raw('D:\S1_data\measurement.dat', Inf, ...
    'LoadIntoMemory', false);

% 第一个成像包对应 StartPulse=1；从成像区间第 1598 个 PRI 读取 30 行
[part, partInfo] = read_s1_raw('', 30, 'StartPulse', 1598, ...
    'CarrierFrequencyHz', 5.405e9, 'Plot', false);
fc = partInfo.radar.fc;                 % Hz
lambda = partInfo.radar.lambda;         % m，299792458/fc
Na_all = partInfo.image.Na_all;         % 始终为全成像区间长度，与本次读 30 行无关
```

- 第二个参数为 **PRI 时隙数**，包含零行；省略或 `Inf` 表示一直读到最后一个成像 ECHO 包的 PRI。有限数量超出末端时也会截到最后一个 ECHO。
- `StartPulse` **从 1 开始，相对全文件首个成像 ECHO 包的 PRI**，默认为 1。它按 PRI 时隙计数，区间内标定和缺包时隙也占一个编号；允许从这些零时隙开始。0、负数、非整数或超出 `Na_all` 会报错。
- `CarrierFrequencyHz` 为正有限标量，默认 `5.405e9`。返回并保存 `radar.fc`（Hz）与 `radar.lambda`（m），不改变原始 I/Q 解码。
- `LoadIntoMemory=true`（默认）返回数值矩阵；`false` 返回 `matfile` 对象，通过 `data.iq(rows, columns)` 读取。
- `Chunk` 参数不再支持，传入会明确报错。所有选定时间范围内的 ECHO chunks 都参与拼接。
- `Plot=true`（默认）显示抽样幅度预览，不改变保存的 I/Q 或零行。
- 可用 `OutputFile` 指定 MAT 路径。每次调用都会实际解码并更新相应输出，输出不得写入原始 SAFE 目录。

2026-10-04 起 `StartPulse` 的基准已改变。换算公式为 `新 StartPulse = 旧 source 零起始偏移 - 首个 ECHO 的 source 零起始偏移 + 1`；本景后者为 528，因此旧 2125 对应新 1598，旧 19800 对应新 19273。新默认文件名使用 `_image_start..._n..._fc...`，避免与旧 `_pri_start...` 输出混淆。

再次使用已生成的结果时无需重新解码。下面的旧文件仍包含 57531 个 source PRI，原有行号、数据和元数据不会因函数升级而自动改变：

```matlab
data = matfile('outputs/LA_S2_HH_full_PRI.mat');
info = load('outputs/LA_S2_HH_full_PRI.mat', 'timeline', 'radar', ...
    'processing_info', 'packet_metadata');
info.mat_file = 'outputs/LA_S2_HH_full_PRI.mat';
part = data.iq(500:2200, :);
% 内存充足时才执行 echo = data.iq;
```

## 本景完整时间轴

产品：`S1A_S2_RAW__0SDH_20150109T014106_20150109T014141_004088_004F1F_2C60.SAFE`，SM 条带成像，S2 波束，HH/HV 双极化。两种极化的 PRI 布局相同。

| 项目 | 数值 |
|---|---:|
| 首个 / 最后成像 PRI Count | 3710 / 60144 |
| 完整成像矩阵 / `image.Na_all` | 56435 × 24424 / 56435 |
| ECHO chunks | 35 |
| 有效回波行 | 55755（每块 1593） |
| 成像区间内标定时隙 | 680 |
| 成像区间内噪声时隙 | 0 |
| 成像区间内 PRI 计数缺口 | 0 |
| 成像区间内总零行 | 680 |
| PRI | 0.000620039222648048 s |
| PRF | 1612.80119624 Hz |
| 距离采样率 | 83410493.8667 Hz |
| 单极化矩阵大小 | 10.27 GiB |

新输出首行就是第一个 ECHO。第一个回波块结束后，第 1594–1613 行是 20 个标定零行，下一回波块从第 1614 行开始。35 块之间共有 680 个标定时隙。整个源文件另有前后非成像区，共 57531 个 PRI（3182–60712），包含 1500 个标定包、16 个噪声包和 260 个无 ISP 时隙；这些全文件统计与 `image.Na_all` 含义不同。

`Na_all` 由**最后成像包的展开 PRI Count − 第一成像包的展开 PRI Count + 1**得到，所以即使成像区间中存在缺包，时隙仍计入。本景成像区间恰好没有无 ISP 的 PRI；代码使用含内部缺包与计数回绕的合成数据验证这一情况。

矩阵行索引对应 PRI Count 的递增位置，而不是简单按包数连续堆叠。缺包会留下零行；PRI Count 的 uint32 正常回绕得到展开；重复、倒退计数或不兼容的采样参数会报错，避免生成错误时间轴。允许 ECHO 的 `Number of Quads`（NQ）与采样窗口长度 SWL 变化，但采样率、子带或 PRI 变化仍会报错。

## 返回参数与文件结构

输出为支持大变量和局部读取的 **MAT v7.3**。Python 分批解码 ECHO 并写入 HDF5，不在内存中创建整景矩阵。每批默认 256 个回波包。写入完成后才替换正式 MAT 文件，并设置 `processing_info.complete=true`；原始 `.dat` 保持只读。

| 变量 | 用途 |
|---|---|
| `iq` | PRI × 距离采样点，complex single |
| `timeline` | 每个矩阵行的 PRI、时间、掩码与对应包/块索引 |
| `radar` | 载频 `fc`、波长 `lambda`、PRI、PRF、采样率及距离向时间参考 |
| `image` | 全成像区间的 `Na_all`、当前行数 `Na` 与首尾成像包的 source PRI 位置 |
| `packet_metadata` | 全文件数据包头，含到矩阵行的映射；枚举使用数字及码表 |
| `chunk_summary` | 全文件 acquisition chunk 概况 |
| `processing_info` | 来源、版本、完成标志和各类时隙数量 |
| `orbit` | 按 GPS 时间去重的 ECEF 位置、速度及时间基准 |
| `attitude` | 按 GPS 时间去重的姿态四元数、本体角速度及时间基准 |

MATLAB 返回的 `info` 包含除 `iq` 外的上述变量，另有 `mat_file` 输出路径。旁边的 `_packets.csv` 与 `_chunks.csv` 便于人工查看。

常用字段：

```matlab
t = info.timeline.time_relative_s;           % 相对全文件第一个 PRI 的秒数
tImage = info.timeline.time_from_image_start_s; % 相对首个成像包的秒数
imagePRI = info.timeline.image_pri_index_one_based; % 相对首个成像包，1 起始
priCount = info.timeline.pri_count_uint32;
isEcho = info.timeline.echo_mask;
isCalibration = info.timeline.calibration_mask;
isNoise = info.timeline.noise_mask;
isMissing = info.timeline.missing_pri_mask;
swst = info.timeline.swst_s;                 % 逐行接收窗口起始时间
numSamples = info.timeline.echo_num_samples; % 原始有效复采样点数：ECHO 为 2*NQ，其他行 0
```

`time_from_selection_start_s` 从当前读取片段的第一个 PRI 起算。`time_relative_s` 和 `global_pri_index_zero_based` 仍保留全文件首包原点，以确保轨道插值时间不变。`packet_index_zero_based`、`global_pri_index_zero_based`、`selected_row_zero_based` 均为零起始索引；缺包没有包索引，使用 -1；所选范围外的包的 `selected_row_zero_based` 也是 -1。`timeline.total_pris` 与 `processing_info.total_pris` 仍是全源文件 PRI 数；全成像区间的长度请使用 `image.Na_all`。

`missing_pri_mask` 表示该 PRI 没有对应 ISP，不直接等同于下传丢包。本景的 Space Packet Count 连续，文件首尾的 PRI 空档符合仪器抑制包生成的情况。

### 轨道与姿态

每次运行 `read_s1_raw` 都从全文件包头的辅助字中提取轨道与姿态，写入 MAT 的 `orbit`、`attitude`，同时返回到 `info.orbit` 和 `info.attitude`。选择 2200 个 PRI 等短片段时仍保留全文件可用的导航记录。

```matlab
[echo, info] = read_s1_raw('', 2200, 'Plot', false);
position = info.orbit.position_ecef_m;           % N x 3，ECEF 位置，米
velocity = info.orbit.velocity_ecef_mps;         % N x 3，ECEF 速度，米/秒
q = info.attitude.quaternion;                   % M x 4，[Q0 Q1 Q2 Q3]
omega = info.attitude.angular_velocity_body_rad_s; % M x 3，本体坐标角速度，rad/s
tOrbit = info.orbit.time_relative_s;
tAttitude = info.attitude.time_relative_s;
```

本景 HH/HV 各有 **37 条轨道记录、37 条姿态记录，时间间隔为 1 秒**。它们按各自时间戳独立排序、去重；行数与回波矩阵不同。位置保留 double 精度，速度、四元数、角速度保留原始 single 精度。

| 字段 | 含义 |
|---|---|
| `available` / `reason` | 是否有可用记录；缺少辅助字段或完整序列时返回 false、原因与固定列数的空数组 |
| `record_count` | 去重后的记录数 |
| `time_gps_s` | GPS 历元 1980-01-06 00:00:00 起算的秒数，**不是 UTC 或 Unix 秒** |
| `time_origin_gps_s` | 全文件第一个包的 `Coarse Time + Fine Time` |
| `time_relative_s` | 相对上述原点的秒数，可与 `info.timeline.time_relative_s` 比较 |
| `orbit.reference_frame` | ECEF；位置和速度均在地心地固坐标系中 |
| `attitude.quaternion_reference_frame` | ECI/J2000；Q0 为标量，Q1–Q3 为向量分量 |
| `attitude.angular_velocity_reference_frame` | spacecraft body，本体坐标系 |

四元数原值、符号和范数均保留，未推定旋转方向；姿态保持原始采样，不插值。斜视角尚未计算。

`read_s1_raw` 还调用 `interpolate_s1_orbit`，把 ECEF 位置和速度计算到所选范围的**每一个 PRI 时隙**，包括标定、噪声和无 ISP 的零行：

```matlab
positionPRI = info.orbit.pulse_position_ecef_m;       % Na x 3，double，m
velocityPRI = info.orbit.pulse_velocity_ecef_mps;     % Na x 3，double，m/s
timePRI = info.orbit.pulse_time_relative_s;          % 与 timeline.time_relative_s 相同
isExtrapolated = info.orbit.pulse_extrapolated_mask; % Na x 1，logical
```

使用位置和速度约束的三次 Hermite 插值，速度取同一插值多项式的导数，在原始导航时间点精确保留原值。短片段采用全文件时间原点，不重复叠加 `StartPulse`。导航缺失时，逐 PRI 数组为 NaN，`pulse_interpolation_available=false` 并记录原因。

这些逐 PRI 字段写入**返回的 `info`**；原始解码 MAT 中仍保存原生导航记录。重新 `load` MAT 后可执行 `info = interpolate_s1_orbit(info)`。SWST 的磁盘输出会连同这些字段保存。

只使用辅助字索引完整为 1–64 且 Space Packet Count 连续的记录，PRI 空档本身不会让序列失效。本景共 923 个序列起点，886 个完整序列，37 个不完整序列跳过；重复时间戳若对应不同数值会报错。相关计数保存在 `complete_subcomm_record_count`、`incomplete_subcomm_record_count` 等字段。

本景导航相对时间约为 −1.156624 至 34.843376 秒；最后一个成像 PRI 超出最后一条导航记录约 **0.475298 秒**。端点外采用最近状态的恒速外推，最多允许一个相邻导航采样间隔，超过会报错。当前全成像区间末端 **767 个 PRI** 由 `pulse_extrapolated_mask` 标记。旧全文件 57531 行时间轴包含尾部非成像区，末端超出约 0.827480 秒、1335 行外推。这里使用的是包头中的星上导航记录，不是外部精密轨道产品。

现有完整 HH 文件 `outputs/LA_S2_HH_full_PRI.mat` 保留原状，已单独生成约 37 KB 的 `outputs/LA_S2_HH_full_PRI_navigation.mat`。已有 `info` 可直接补入，无需重新读取回波：

```matlab
nav = load('outputs/LA_S2_HH_full_PRI_navigation.mat');
info.orbit = nav.orbit;
info.attitude = nav.attitude;
info = interpolate_s1_orbit(info);
```

新运行生成的 MAT 已直接包含两项变量；旧文件不会自动升级。时间、单位和参考系依据 [Sentinel-1 SAR Space Packet Protocol §3.2.3](https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf#page=23)。

### 距离窗口的变化

`read_s1_raw` 以**全文件所有 ECHO 包中最大的 `2*NQ`** 作为输出矩阵列数，与 `StartPulse` 和本次读取长度无关。每个 ECHO 包解码出的原始复采样点从该行第 1 列写入；较短回波仅在末尾补零，不缩放、不重采样。NQ/SWL 的变化不会改变完整 PRI 时间轴，标定、噪声和缺包行仍保持全零。

`info.timeline.echo_num_samples` 与所选矩阵行一一对应，记录各 ECHO 包的原始有效复采样点数 `2*NQ`；非 ECHO 和无包时隙记为 0。它表示解码样本长度，不按样本是否非零来计数，也不包含末尾补零。不同 SWST 的回波仍各自以原接收窗口起点为第一个样本，统一距离起点需另外调用 `echo_SWST`。

本景 ECHO 的 SWST 有两种值，差值对应 **100 个距离采样点**，首次变化位于 PRI 23066（新成像编号第 19357 行，旧全文件 MATLAB 行 19885）。`read_s1_raw` 返回原始接收窗口内的采样排列，`radar.shared_absolute_range_origin=false`，单值 `radar.swst_s=NaN`，逐行值见 `timeline.swst_s`。没有 ISP 的行原始 SWST 为 NaN。

### IQ 校正

`IQ_correction.m` 按参考 C++ MEX 源码的 `iqBiasCorrection` 公式，先逐脉冲去除 I/Q 均值、将 Q 乘以 `std(I)/std(Q)`，再按本次提取数据中 `swath_number` 的连续变化分段。每段独立对有效脉冲相关系数作 Fisher z 等权平均，估计该段统一的非正交角并修正 Q。它仅支持内存中的 `single/double` 回波矩阵，保持原尺寸和精度。

```matlab
[echo, info] = read_s1_raw(inputFile, Na, 'StartPulse', startPulse, 'Plot', false);
[echo, info] = IQ_correction(echo, info);  % 在 SWST 对齐之前执行
[echo, info] = echo_SWST(echo, info);
```

子带号优先取 `info.timeline.swath_number`，没有该字段时，从 `packet_metadata.values.Swath_Number` 按 `selected_row_zero_based` 映射到当前矩阵行。`1→2→1` 形成三个独立段，不合并两段编号 1；标定和噪声包的编号变化也用于确定边界。缺包不单独开启新段，其段号沿用前一段；选区开头无包时段号为 0，直到遇到首个已知编号。

只有 `echo_mask=true` 行的 `1:echo_num_samples` 参与统计和校正；标定、噪声、缺包零行以及短回波尾部补零保持不变。零方差或近乎完全相关的行仅去均值，不参与角度估计；整段无可估计回波时该段角度记为 0。`info.iq_correction` 保存逐行 `dc_i`、`dc_q`、`gain_q`、`estimated_mask`、`skipped_mask`，以及下列分段结果：

```matlab
info.iq_correction.segment_id              % 每行所属段号；开头未知段为 0
info.iq_correction.segment_start_row       % 各段在当前 echo 中的起始行（从 1 开始）
info.iq_correction.segment_end_row         % 各段结束行
info.iq_correction.segment_swath_number    % 各段子带号
info.iq_correction.theta_rad               % 各段非正交角，列向量，单位 rad
info.iq_correction.segment_estimated_count % 各段参与角度估计的脉冲数
```

该函数不写回原始 MAT 文件，`Read_Data_RD.m` 已按上述顺序调用。

### 按 SWST 对齐距离采样

```matlab
% 内存中的 echo：返回新的 complex single/double 矩阵。
[echoAligned, infoAligned] = echo_SWST(echo, info);

% 整景推荐磁盘模式：分批处理，返回新的 matfile（变量仍名为 iq）。
[data, info] = read_s1_raw('', Inf, 'LoadIntoMemory', false, 'Plot', false);
[dataAligned, infoAligned] = echo_SWST(data, info);
block = dataAligned.iq(1:256, :);
disp(infoAligned.mat_file);  % 新的 *_SWST_*.mat 路径
```

已生成的完整 HH 文件可以使用上面的 `matfile`、元数据 `load` 和导航补入示例，然后直接执行 `echo_SWST(data, info)`，无需重新解码。磁盘模式每批 64 行，原始输入文件不变；输出完成后才发布新文件。输入的行与 `info.timeline` 必须一一对应。对已对齐的数据和对应 `info` 再次调用将直接返回，避免重复时移。

只用 `echo_mask=true` 的有效回波行计算最小、最大 SWST，标定和无包行不参与窗口范围估计，始终保持全零。设 `Fs` 为距离采样率，原长度为 `Nr0`：

```matlab
delta = (SWST - minSWST)*Fs;                 % 逐行右移的采样点数
Nr = Nr0 + ceil((maxSWST-minSWST)*Fs);
Nr = Nr + mod(Nr, 2);                       % 补成偶数
k = [0:Nr/2-1, -Nr/2:-1];
% 每行末尾先补零，再乘线性相位，支持分数采样点时移：
y = ifft(fft(x, Nr).*exp(-1i*2*pi*k*delta/Nr));
```

整数附近仅由浮点误差引起的偏差会先归整，避免多补一点。分数时移采用周期带限插值，可能产生 sinc 尾部；标定等零行仍精确保留为零。本景整景或跨 SWST 变化的片段从 **24424 列扩为 24524 列**；若所选片段内 SWST 相同，则只检查原长度的奇偶性。

对齐后 `timeline.swst_s` 与 `radar.swst_s` 为统一的 `minSWST`，`timeline.range_window_offset_samples` 全零。原始逐行 SWST 保存在 `swst_alignment.original_swst_s`，实际时移见 `swst_alignment.shift_samples`，补零数量见 `swst_alignment.padding_samples`。`packet_metadata` 保持原始包头数值。

`timeline.echo_num_samples` 在 SWST 对齐后仍保留原始样本长度。由于各行起点已平移，它不再等同于该行有效样本的终止列号；分数时移还可能产生插值尾部，因此不能仅据此字段截断对齐后的回波。

统一距离轴和中心采样点：

```matlab
range_m = infoAligned.radar.slant_range_m;                 % Nr x 1，m
centerIndex = infoAligned.radar.center_sample_index;      % Nr/2+1
centerRange_m = infoAligned.radar.center_slant_range_m;   % range_m(centerIndex)
```

`radar.sampling_window_time_s` 以 `minSWST` 为起点。转换为 RAW 回波斜距还计入包头 `Rank*PRI` 和协议规定的抽取滤波瞬态丢弃时间：

`R(n) = c/2 * (Rank*PRI + minSWST + 320/(8*f_ref) + n/Fs)`，`n=0,...,Nr-1`，`f_ref=37.53472224 MHz`。

各量分别存入 `radar.slant_range_time_s`、`decimation_transient_delay_s` 和 `range_sample_spacing_m`。本景两种 SWST 都包含时，距离间隔约 **1.797090774 m**，首点约 **754754.163957 m**，第 **12263** 点（`Nr/2+1`）的中心距离约 **776790.091031 m**。该值是未聚焦 RAW 的采样时序斜距，尚未应用 SWST bias 标定或距离压缩后的参考点修正；`swst_bias_applied=false`，`swst_bias_s=NaN`。依据 [Sentinel-1 SAR Space Packet Protocol §3.2.5.9–11](https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf)。

### 漂移与内部时延补偿

`Read_Data_RD.m` 已按以下顺序调用，回波为内存中的 `Na×Nr` 矩阵：

```matlab
[echo, info] = IQ_correction(echo, info);
[echo, info] = echo_SWST(echo, info);
[echo, info] = Shift_Delay_correction(echo, info, inputFile);
```

第三个参数可省略，默认读取 `info.processing_info.source_file`。函数在 RAW 产品的同级 `*SLC*` 产品目录下查找 `*.SAFE/annotation` 主 XML，按卫星、子带、极化、开始时刻、轨道号和 datatake 匹配，允许 RAW/SLC 结束时刻和校验码不同。当前支持单一子带；没有对应 XML 或匹配不唯一时直接报错。

这两项未找到 MEX 实现，采用参考 MATLAB 的 `instrument_drift_compensation.m`、`internal_calibration.m` 和 `range_compression.m` 公式：

- 漂移：将有效 replica 的 `pgProductAmplitude * exp(1i*pgProductPhase)` 按脉冲时刻做复数 `spline` 插值，再乘到对应回波行。按参考源码允许外推，并在 `extrapolated_mask` 中标记。RAW 的 GPS 时刻先通过 MATLAB 闰秒表转为 UTC，与 XML 时间对齐。
- 内部时延：取 `referenceReplica.timeDelay`（秒），沿距离维乘 `exp(-1i*2*pi*f*delay)`。正值向右移动；不再加入已由 `echo_SWST` 处理的 SWST 项。当前 S1/VH XML 为 `4.330593e-7 s`，约 `43.3460` 个采样点。

内部时移先补零做 FFT，再截回原来的 `Nr` 列，因此矩阵尺寸、距离网格和中心距离字段不变；移出窗口的信号被截去，分数时移可能出现插值尾部。标定和缺包行始终为零，原始有效长度保留。补偿参数写入 `info.shift_delay_correction`，其中 `pg_product` 为逐行复数系数，`internal_time_delay_s` 为秒，`internal_delay_samples` 为采样点。另有 `info.radar.internal_delay_applied=true`；后续距离压缩不要再次补偿这项内部时延。SWST bias 标记保持独立。函数不写回原始 MAT。

### 多普勒与几何参数

`Read_Data_RD.m` 在 `Shift_Delay_correction` 后增加 `info = doppler_geo_parameter(info);`。函数只读取上一步已匹配的 SLC XML 并更新 `info`，不修改 `echo`。C++/MEX 中未找到对应实现，主要依据 MATLAB 的 `parse_xmlpara`、`get_doppler_parameter`、`get_target_pos` 和 `get_equivalent_vel`。

| 字段 | 内容与单位 |
| --- | --- |
| `info.doppler.centroid_hz` | 当前选段代表性的绝对多普勒中心，Nr×1，Hz |
| `info.doppler.geometry_centroid_hz / data_centroid_hz` | XML 几何/数据多项式各自的结果，Nr×1，Hz |
| `info.doppler.azimuth_fm_rate_hz_s` | 方位调频率，Nr×1，Hz/s |
| `info.doppler.centroid_rate_hz_s` | 源码按方位波束扫描率计算的多普勒中心变化率，Hz/s |
| `info.doppler.processing_bandwidth_hz` | XML 方位处理带宽，Hz，并非独立估计的瞬时带宽 |
| `info.xml_orbit` | XML 原始 UTC/ECEF 轨道与插值到每个 PRI 的位置、速度，分别为 m、m/s |
| `info.geometry.xml_reference_range_m` | XML `referenceRange`：距离扩散损失归一化的参考距离，m |
| `info.geometry.center_slant_range_m` | 当前 RAW 矩阵中心采样点的斜距，m |
| `info.geometry.target_position_ecef_m / target_lla_deg_m` | 中心目标的 ECEF 位置和 `[纬度deg, 经度deg, 高程m]` |
| `info.geometry.target_sample_position_ecef_m` | 按距离门采样计算的目标位置，Ns×3；列号见 `sample_range_indices` |
| `info.geometry.equivalent_velocity_mps` | 所有距离门的等效速度，Nr×1，m/s |
| `info.geometry.satellite_speed_mps` | 每个 PRI 的卫星 ECEF 速率，Na×1，m/s |

多普勒记录选择沿用源码：分别选取与所选完整 PRI 时间轴总时间差最小的一条 DC/FM 记录，用完整双程距离时间 `info.radar.slant_range_time_s - t0` 求多项式。系数、记录索引/时间、数据质量标志、原始 XML DC/FM 记录均保存在 `info.doppler`。这是一组代表当前选段的距离向参数，不是逐脉冲或逐 TOPS burst 的不同多项式。

修正参考源码重复解模糊的问题：SLC 的 `dataDcPolynomial` 已表示绝对多普勒中心，因此直接求值；数据质量超阈值时使用几何多项式，并记录 `centroid_source`。源码的 `fix(geometryDC/PRF)*PRF + dataDC` 仍保存在 `reference_formula_centroid_hz` 供比较。当前场景两者一致。依据 [ESA Level-1 算法 §5.4–5.5](https://sentinels.copernicus.eu/documents/247904/1877131/Sentinel-1-Level-1-Detailed-Algorithm-Definition) 与 [TOPS SLC Deramping Eq.13](https://sentiwiki.copernicus.eu/__attachments/1673968/COPE-GSEG-EOPG-TN-14-0025%20-%20Sentinel-1%20TOPS%20SLC%20Deramping%202017%20-%201.3.pdf)。

几何计算优先使用现有 RAW 逐脉冲轨道；不可用时采用 XML 轨道，实际来源见 `info.geometry.orbit_source`。XML 六个轨道分量按源码独立 spline 插值，外推由 `info.xml_orbit.extrapolated_mask` 标记，不覆盖 `info.orbit`。中心行取 `floor(Na/2)+1`，中心列与现有距离轴一致；每隔 30 个距离门求解零多普勒、零高程、右视 WGS84 目标，再根据所选脉冲的距离历程拟合等效速度，沿距离 spline 插值。额外纳入中心列和最后一列。该位置是源码的椭球模型近似，不使用 DEM 或实测多普勒定位；短时间段的等效速度更易受轨道精度影响。

其他参数包括最近斜距/对应离散时刻、定位约束残差、方位转向率、XML 载频、SLC 首采样距离时间等。XML 载频和距离参数不替换用户设置的载频及 RAW 距离轴。源码固定天线尺寸和波束宽度放在 `geometry.reference_defaults`；`beam_command_squint_rad` 仅为包头波束地址换算的指令角，不能当作从多普勒反算的斜视角。`xml_zero_doppler_minus_acquisition_s` 只保留 XML 元数据，不额外平移上述源码模型的时刻。

## 环境与验证

本目录 `.venv` 为 Python 3.12 环境，依赖固定在 `requirements.txt`。MATLAB 使用外部 Python 进程，不依赖 `pyenv`。在另一台 Windows 电脑配置：

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

可通过 `'Python', 'D:\envs\s1\python.exe'` 指定其他兼容环境。

验证入口：

- `test_decode_s1_raw.py`：合成数据检查 PRI 缺口、计数回绕、跨 chunk 放置、零行、SWST 和不完整输出保护。
- `test_decode_s1_navigation.py`：用实际大端字节布局构造辅助字，检查轨姿解码、独立时间排序去重、截断/断包、计数回绕、异常输入和 MAT 字段形状。
- `validate_read_s1_raw.m`：实际 HH/HV 各读取从首 ECHO 开始的 2200 个 PRI；检查内部标定零行、StartPulse=1/2 的参考回波、尾部 Inf/超长请求截断、自定义载频与 Na_all 的 MAT 持久化；另验证标定间隔截取、逐 PRI 轨道，并在 MATLAB 独立重组辅助字比对导航数值。
- `validate_full_pri_matrix.m`：只读既有完整 HH MAT，检查首尾 ECHO 界定的 56435 行成像区间、零行及原参考回波，不重新解码或覆盖旧文件。
- `validate_s1_orbit_interpolation.m`：解析三次轨迹、原始节点、选段时间基准、外推边界和异常输入；只读取整景元数据，验证 57531 行逐 PRI ECEF 轨道。
- `validate_echo_SWST.m`：整数/分数时移、偶数长度、零时隙、距离轴、重复调用及新 MAT 输出；实际读取 SWST 变化前后的 256 个 PRI 并对照直接整数移位。
- `validate_s1_variable_echo.m`：实际读取 2018 S1/VH 数据的 1024 脉冲、长度变化交界和仅短回波选区，检查尾补零、原始长度持久化、磁盘读取及 SWST 兼容；同时检查原 S2 数据的等长回波和标定零行。
- `validate_IQ_correction.m`：已知 I/Q 失衡恢复、single/double、变长零尾与零时隙、逐段等脉冲权重、子带返场/标定边界/缺包、包头选区映射、退化样本及输入检查；只读已有实际回波验证 IQ 校正后再进行 SWST 对齐。
- `validate_Shift_Delay_correction.m`：正负时移与边界、复数 PG 插值、零行、单双精度、GPS/UTC 闰秒转换和真实 32 脉冲的 IQ→SWST→漂移/内部时延流程。
- `validate_doppler_geo_parameter.m`：已知多项式、绝对 DC/质量标志、XML 轨道、解析右视目标、已知等效速度、奇偶尺寸，以及真实 32/1024 脉冲完整流程。

旧的 `LA_S2_HH_first256.mat` / `LA_S2_HV_first256.mat` 是首个 ECHO 块的 256 行参考结果，不包含文件开头非回波时隙，保留用于回归对照。

2026-09-08 已在 MATLAB R2024a Update 8 实际通过上述 HH/HV 跨块、标定区间截取和完整 HH 矩阵验证；全部 57531 行扫描通过，55755 行回波均有有效值，1776 行非回波/缺包位置全部为零。8 个合成测试也通过。实际日志为 `outputs/matlab_pri_subset_validation.log` 和 `outputs/matlab_full_pri_validation.log`。

轨道与姿态加入后，16 项 Python 合成测试通过；MATLAB R2024a 实际 HH/HV 调用均返回 37 条轨道与 37 条姿态记录，坐标、时间、四元数和角速度通过独立辅助字重组比对。独立导航文件与新调用结果一致，原回波参考片段保持逐元素一致。日志为 `outputs/matlab_navigation_validation.log`。

逐 PRI 轨道集成后，实际 HH/HV 2200 行和 30 行选段检查通过，选段 ECEF 与全景对应行完全一致，日志为 `outputs/matlab_orbit_reader_validation.log`。SWST 验证通过：实际 256 行含 236 行回波、20 行标定，输出 24524 列；FFT 时移与直接右移 100 点的相对误差约 `3.99e-7`，内存与磁盘输出逐元素一致，输入保持原值。日志为 `outputs/matlab_swst_validation.log`。

2026-10-04 新接口通过 14 项读取器合成测试与 8 项导航测试。MATLAB R2024a 实际 HH/HV 首 2200 行、成像起点 1/2、6 GHz 载频、自定义选段及末尾裁切均通过，`Na_all` 始终为 56435，source/image 时间轴与新字段持久化正确。新 StartPulse=19273 的 SWST 变化片段仍输出 256×24524，误差约 `3.99e-7`。既有整景 MAT 保持只读，修正验证脚本的磁盘索引方式后，代表段中的 624 行回波和 40 行零行也检查通过。此次没有重新解码整景。

2026-10-04 变长回波支持通过 17 项读取器测试与 8 项导航测试。实际 S1/VH 新景读取 `StartPulse=1024, Na=1024` 得到 `1024×25780`；跨第 3743 个成像脉冲的 32 行包含 16 行 25780 点和 16 行 25764 点，短行尾部 16 点为零，有效 I/Q 与直接解码结果逐元素一致。MATLAB 验证原始长度字段为逐行 `uint32` 列向量，跨界 SWST 对齐输出 25912 列且保留原始长度；仅短回波的磁盘读取仍为 25780 列。原 S2 28 行验证包含 20 行标定零行与 8 行 24424 点回波，SWST 合成和磁盘读写测试也通过。日志为 `outputs/matlab_variable_echo_validation.log`。

2026-10-05 漂移与内部时延补偿通过 MATLAB R2024a 验证：正负整点/分数时移、独立解析核、single/double、零行、复数 PG 三次插值与外推、有效标志过滤，以及 2015/2018 的 GPS→UTC 时间转换均通过。真实 S1/VH 32 脉冲完成 IQ→SWST→漂移/内部时延处理，保持 `32×25912`，内部时移 `43.3460281` 点。日志为 `outputs/matlab_shift_delay_validation.log`。

2026-10-05 多普勒与几何参数通过 MATLAB R2024a 验证：合成多项式、绝对 DC 防重复解模糊、质量标志回退、XML 轨道及 RAW 优先、奇偶尺寸、解析右视椭球目标与已知 7500 m/s 等效速度均通过。真实 32/1024 脉冲完整流程通过，中心等效速度分别约 `7195.529015 / 7196.539187 m/s`。日志为 `outputs/matlab_doppler_geo_validation.log`。

RAW I/Q 尚未进行距离压缩、方位聚焦、整幅图像地理编码或完整辐射定标。acquisition chunk 是解码库的包分组，不等同于 TOPS burst。

相关资料：[sentinel1decoder](https://github.com/Rich-Hall/sentinel1decoder)、[ESA Sentinel-1 产品说明](https://sentiwiki.copernicus.eu/web/s1-products)、[MATLAB 局部加载 MAT 文件](https://www.mathworks.com/help/matlab/import_export/load-parts-of-variables-from-mat-files.html)。核心解码由第三方开源库完成，授权条款见其仓库。
