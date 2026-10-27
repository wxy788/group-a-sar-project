function [echo, info] = IQ_correction(echo, info)
%IQ_CORRECTION 逐行去均值、平衡 I/Q 增益，再按连续子带段校正非正交误差。
% 输入 echo：行是脉冲，列是距离采样点；info 来自 read_s1_raw。
% 顺序：read_s1_raw -> IQ_correction -> echo_SWST。
% 只处理原始有效样本，零行和尾部补零保持不变。
% info.iq_correction 保存校正参数；theta_rad 为逐段角度，单位 rad。

validateattributes(echo, {'single','double'}, {'2d','nonempty','nonsparse'});
assert(isstruct(info) && isfield(info, 'timeline') && ...
    all(isfield(info.timeline, {'echo_mask', 'echo_num_samples'})), ...
    'IQ_correction:MissingInfo', '请使用新版 read_s1_raw 返回的 info。');
assert(~isfield(info, 'swst_alignment') || ~info.swst_alignment.applied, ...
    'IQ_correction:AlreadyAligned', '请在 echo_SWST 之前调用 IQ_correction。');
assert(~isfield(info, 'iq_correction') || ~info.iq_correction.applied, ...
    'IQ_correction:AlreadyCorrected', '该回波已进行 IQ 校正，请勿重复调用。');

%% 1. 读取每行信息，按 swath_number 的连续变化分段
[Na, Nr] = size(echo);
echoMask = info.timeline.echo_mask(:);
numSamples = double(info.timeline.echo_num_samples(:));
assert(numel(echoMask) == Na && all(ismember(echoMask, [0, 1])) && ...
    numel(numSamples) == Na && all(isfinite(numSamples)) && ...
    all(numSamples == fix(numSamples) & numSamples >= 0 & numSamples <= Nr), ...
    'IQ_correction:InvalidLengths', '逐行回波标记、有效长度必须与 echo 尺寸一致。');
echoMask = logical(echoMask);
assert(all(numSamples(echoMask) > 0) && all(numSamples(~echoMask) == 0), ...
    'IQ_correction:InvalidLengths', '回波行长度应大于零，非回波行长度应为零。');
rows = find(echoMask).';
swath = getSwathNumber(info, Na);
assert(all(isfinite(swath(echoMask))), 'IQ_correction:MissingSwath', ...
    '每个回波行都需要对应的 swath_number。');

% 例如 1→2→1 分成三段；标定包也确定边界，缺包沿用前段编号。
knownRows = find(isfinite(swath));
segmentId = zeros(Na, 1);
startRows = zeros(0, 1); endRows = zeros(0, 1);
if ~isempty(knownRows)
    startRows = knownRows([true; diff(swath(knownRows)) ~= 0]);
    endRows = [startRows(2:end)-1; Na];
    segmentId(startRows) = 1;
    segmentId = cumsum(segmentId);
end

%% 2. 逐行估计均值、Q 通道增益比、I/Q 相关系数
dc = complex(nan(Na, 1), nan(Na, 1));      % 实部是 I 均值，虚部是 Q 均值
gainQ = nan(Na, 1); z = nan(Na, 1);
for k = rows
    x = double(echo(k, 1:numSamples(k)));   % 仅取有效样本，用 double 统计
    assert(all(isfinite(x)), 'IQ_correction:NonfiniteEcho', ...
        '第 %d 行有效回波含 NaN 或 Inf。', k);
    dc(k) = mean(x);
    x = x - dc(k);
    I = real(x); Q = imag(x);
    sigmaI = sqrt(mean(I.^2));             % 总体标准差，分母为有效样点数
    sigmaQ = sqrt(mean(Q.^2));
    correlation = mean(I.*Q) / (sigmaI*sigmaQ);
    gainQ(k) = 1;
    % 零方差或完全相关时无法估计相位，这些行仅去均值。
    if isfinite(correlation) && abs(correlation) < 1 - 1e-12
        gainQ(k) = sigmaI / sigmaQ;
        z(k) = atanh(correlation);         % Fisher 变换，供下一步按段平均
    end
end

%% 3. 每段估计一个非正交角，再用所在段的角度校正各行
estimated = isfinite(z);
theta = zeros(numel(startRows), 1);
estimatedCount = zeros(numel(startRows), 1);
for s = 1:numel(startRows)
    segmentRows = startRows(s):endRows(s);
    usedRows = segmentRows(estimated(segmentRows));
    estimatedCount(s) = numel(usedRows);
    if ~isempty(usedRows)
        theta(s) = asin(tanh(mean(z(usedRows))));
    end                                   % 无可用回波的段保留 theta=0
end

for k = rows
    cols = 1:numSamples(k);
    x = double(echo(k, cols)) - dc(k);
    I = real(x); Q = imag(x);
    if estimated(k)
        angle = theta(segmentId(k));
        Q = gainQ(k)*Q;
        Q = Q/cos(angle) - I*tan(angle);
    end
    echo(k, cols) = complex(I, Q);         % 只改有效区间，保留输入精度
end

info.iq_correction = struct('applied', true, 'dc_i', real(dc), 'dc_q', imag(dc), ...
    'gain_q', gainQ, 'theta_rad', theta, 'estimated_mask', estimated, ...
    'skipped_mask', echoMask & ~estimated, 'segment_id', segmentId, ...
    'segment_start_row', startRows, 'segment_end_row', endRows, ...
    'segment_swath_number', swath(startRows), 'segment_estimated_count', estimatedCount);
end

function swath = getSwathNumber(info, Na)
% 包头表包含整个文件，必须用选区行号映射，不能直接取前 Na 个包。
if isfield(info.timeline, 'swath_number')
    swath = double(info.timeline.swath_number(:));
else
    assert(isfield(info, 'packet_metadata') && isfield(info.packet_metadata, 'values') && ...
        all(isfield(info.packet_metadata.values, {'Swath_Number', 'selected_row_zero_based'})), ...
        'IQ_correction:MissingSwath', 'info 中缺少逐行子带号或包头到选区的映射。');
    packets = info.packet_metadata.values;
    row = double(packets.selected_row_zero_based(:));
    values = double(packets.Swath_Number(:));
    assert(numel(row) == numel(values) && isreal(row) && ...
        all(isfinite(row) & row == fix(row) & row >= -1 & row < Na), ...
        'IQ_correction:InvalidSwath', '包头选区行号无效。');
    selected = row >= 0;                   % -1 表示本次未提取的包
    assert(numel(unique(row(selected))) == nnz(selected), ...
        'IQ_correction:InvalidSwath', '多个包映射到了同一个回波矩阵行。');
    swath = nan(Na, 1);                    % 无包的 PRI 没有子带号
    swath(row(selected)+1) = values(selected);
end
known = ~isnan(swath);
assert(numel(swath) == Na && isreal(swath) && ...
    all(isfinite(swath(known)) & swath(known) == fix(swath(known)) & ...
        swath(known) >= 0 & swath(known) <= 255), ...
    'IQ_correction:InvalidSwath', '子带号应逐行对应，取 0~255 整数；缺包可为 NaN。');
end
