function [echo, info] = Shift_Delay_correction(echo, info, inputFile)
%SHIFT_DELAY_CORRECTION 在 echo_SWST 后补偿设备漂移和内部时延。
% [echo, info] = Shift_Delay_correction(echo, info, inputFile)
% inputFile 为 RAW .dat 路径，省略时使用 info.processing_info.source_file。
% 自动查找同级 SLC 产品的主 annotation XML，适用于单一子带的回波矩阵。
% 漂移：逐行乘插值后的 PG；时延：距离 FFT 乘 exp(-1i*2*pi*f*delay)。
% 输入输出均为 Na×Nr，保持 single/double 类型，标定和缺包行保持全零。
% 时移先补零，最后截回原 Nr 列；移出窗口的信号被截去，分数时移有插值尾部。
% 参考 MATLAB 源码：instrument_drift_compensation / internal_calibration /
% range_compression。SWST 已完整对齐，这里不再重复补偿 SWST。

if nargin < 3, inputFile = info.processing_info.source_file; end
validateattributes(echo, {'single','double'}, {'2d','nonempty','nonsparse','finite'});
assert(isfield(info, 'swst_alignment') && info.swst_alignment.applied, ...
    'Shift_Delay_correction:NotAligned', '请先调用 echo_SWST。');
assert(~isfield(info, 'shift_delay_correction') || ~info.shift_delay_correction.applied, ...
    'Shift_Delay_correction:AlreadyCorrected', '已完成漂移和内部时延补偿，请勿重复调用。');
[Na, Nr] = size(echo);
echoMask = info.timeline.echo_mask(:);
relativeTime = double(info.timeline.time_relative_s(:));
assert(numel(echoMask) == Na && all(ismember(echoMask, [0,1])) && ...
    numel(relativeTime) == Na && all(isfinite(relativeTime)), ...
    'Shift_Delay_correction:RowMismatch', 'echo 与 info 的逐脉冲信息必须对应。');
echoMask = logical(echoMask);
assert(all(echo(~echoMask,:) == 0, 'all'), ...
    'Shift_Delay_correction:NonzeroInactiveRow', '标定及缺包行应为零。');
Fs = double(info.radar.range_sample_rate_hz);
validateattributes(Fs, {'double'}, {'scalar','real','finite','positive'});

%% 1. 找到相同观测、相同极化的 SLC XML
[xmlFile, swath] = findAnnotation(inputFile);
xml = readstruct(xmlFile);
replicaInfo = xml.generalAnnotation.replicaInformationList.replicaInformation;
replicaInfo = replicaInfo(strcmpi(string({replicaInfo.swath}), swath));
assert(isscalar(replicaInfo), 'Shift_Delay_correction:SwathMismatch', ...
    'XML 必须包含当前子带唯一的 replicaInformation。');
delay = double(replicaInfo.referenceReplica.timeDelay);  % 秒，不叠加各 replica 的时延
validateattributes(delay, {'double'}, {'scalar','real','finite'});
assert(abs(delay*Fs) < Nr, 'Shift_Delay_correction:DelayTooLarge', ...
    '内部时延已超出整个距离窗口，请核对 XML 与采样率。');
replica = replicaInfo.replicaList.replica;
valid = strcmpi(string([replica.relativePgProductValidFlag]), 'true') & ...
    strcmpi(string([replica.absolutePgProductValidFlag]), 'true');
replica = replica(valid);
assert(~isempty(replica), 'Shift_Delay_correction:NoValidPG', 'XML 中没有有效的 PG 系数。');
replicaTime = [replica.azimuthTime].';
if ~isdatetime(replicaTime)
    replicaTime = datetime(replicaTime, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");
end
replicaTime.TimeZone = 'UTC';
pg = double([replica.pgProductAmplitude].') .* exp(1i*double([replica.pgProductPhase].'));
[replicaTime, order] = sort(replicaTime);
pg = pg(order);
assert(all(isfinite(pg)) && all(~isnat(replicaTime)) && ...
    all(diff(replicaTime) > seconds(0)), 'Shift_Delay_correction:InvalidPG', ...
    'PG 系数和时刻必须有效，时刻不能重复。');

%% 2. 完整 PRI 时间轴由 GPS 转成 UTC，再插值复数 PG
% 相对时间从全文件第一个包起算，已包含 StartPulse 偏移，不要再加一次。
if isfield(info, 'orbit') && isfield(info.orbit, 'time_origin_gps_s')
    origin = double(info.orbit.time_origin_gps_s);
else
    packet = info.packet_metadata.values;
    origin = double(packet.Coarse_Time(1)) + double(packet.Fine_Time(1));
end
validateattributes(origin, {'double'}, {'scalar','real','finite'});
pulseTime = datetime(1980,1,6, 'TimeZone','UTCLeapSeconds') + seconds(origin + relativeTime);
pulseTime.TimeZone = 'UTC';                            % MATLAB 自动处理闰秒
sourceTime = seconds(replicaTime - replicaTime(1));
targetTime = seconds(pulseTime - replicaTime(1));
if isscalar(pg)
    pgPulse = repmat(pg, Na, 1);
else
    pgPulse = interp1(sourceTime, pg, targetTime, 'spline', 'extrap');
end
assert(all(isfinite(pgPulse)), 'Shift_Delay_correction:InvalidPG', 'PG 插值结果无效。');

%% 3. 按行乘 PG，在距离频域补偿内部时延
% 与参考源码同号：delay>0 向右移。只移动数据，不再平移距离轴。
Nfft = 2^nextpow2(Nr + ceil(abs(delay*Fs)));
frequency = ifftshift(-floor(Nfft/2):ceil(Nfft/2)-1) * (Fs/Nfft);
phase = cast(exp(-1i*2*pi*frequency*delay), 'like', echo);
rows = find(echoMask).';
for first = 1:32:numel(rows)
    selected = rows(first:min(first+31, numel(rows)));
    block = echo(selected,:) .* cast(pgPulse(selected), 'like', echo);
    if delay ~= 0
        block = ifft(fft(block, Nfft, 2).*phase, [], 2);
    end
    echo(selected,:) = block(:,1:Nr);
end

%% 4. 保存实际采用的参数，避免后续距离压缩重复补偿
info.shift_delay_correction = struct('applied', true, 'xml_file', xmlFile, ...
    'pg_product', pgPulse, 'internal_time_delay_s', delay, ...
    'internal_delay_samples', delay*Fs, 'fft_length', Nfft, ...
    'pulse_time_posix_s', posixtime(pulseTime), ...
    'replica_time_posix_s', posixtime(replicaTime), 'replica_pg_product', pg, ...
    'extrapolated_mask', targetTime < sourceTime(1) | targetTime > sourceTime(end), ...
    'row_applied_mask', echoMask, 'pg_source', char(replicaInfo.referenceReplica.pgSource), ...
    'method', 'Complex PG spline/extrap; range FFT exp(-1i*2*pi*f*delay)', ...
    'boundary_note', 'Zero-padded FFT, cropped to original Nr; fractional shifts have interpolation tails.');
info.radar.internal_time_delay_s = delay;
info.radar.internal_delay_applied = true;
% 距离网格及中心采样点不变；SWST bias 与内部时延是不同字段，不改其标记。
info.mat_file = '';                                    % 当前结果只在内存中
end

function [xmlFile, swath] = findAnnotation(inputFile)
% RAW 和 SLC 的结束时刻、校验码可以不同，按开始时刻/轨道/数据条带/极化匹配。
inputFile = char(inputFile);
[safeDir, rawName] = fileparts(inputFile(:).'); % Python 写入的 MAT 字符串可能是列向量
parts = regexp(lower(rawName), ...
    '^(s1[a-z])-([a-z0-9]+)-raw-[a-z]-([hv]{2})-(\d{8}t\d{6})-\d{8}t\d{6}-(\d{6})-([0-9a-f]{6})$', ...
    'tokens', 'once');
assert(~isempty(parts), 'Shift_Delay_correction:InputFile', 'inputFile 应为 Sentinel-1 RAW .dat 路径。');
dataDir = fileparts(fileparts(safeDir));
pattern = sprintf('%s-%s-slc-%s-%s-*-%s-%s-*.xml', parts{:});
files = dir(fullfile(dataDir, '*SLC*', '*.SAFE', 'annotation', pattern));
assert(isscalar(files), 'Shift_Delay_correction:AnnotationNotUnique', ...
    '应找到唯一对应的 SLC annotation XML，实际找到 %d 个。查找目录：%s', numel(files), dataDir);
xmlFile = fullfile(files.folder, files.name);
swath = upper(parts{2});
end
