function read_level0(varargin)
%READ_LEVEL0  Sentinel-1 Level-0 RAW 解码与复数读取封装（A组·第2周交付）
%
% 功能：调用助教工具包 sentinel_l0_toolkit/read_s1_raw 解码 RAW 回波，
%       可选做 IQ 校正与 SWST 距离窗对齐，保存复数块、参数(JSON) 与 3 个 ROI。
%
% 用法：
%   read_level0('RawRoot', rawRoot)                                % 默认读一段 4000 行块
%   read_level0('RawRoot', rawRoot, 'Full', true)                  % 整景解码（磁盘模式，~10GB）
%   read_level0('RawRoot', rawRoot, 'Polarization','VV', ...
%               'StartPulse',1598, 'NumPulses',4000, ...
%               'OutDir','D:/study/SAR/data/processed/level0', ...
%               'ROIs', {[r1 c1 h w], [r2 c2 h w], [r3 c3 h w]})
%
% 输入参数：
%   RawRoot        RAW 产品根目录（含 *SAFE 的文件夹），必填
%   Polarization   'VV' 或 'VH'（本产品 SDV 双极化；非 HH），默认 'VV'
%   StartPulse     相对首个 ECHO 的 1 起始 PRI，默认 1
%   NumPulses      内存块行数（含零行），默认 4000；整景请用 Full=true
%   Full           true=整景磁盘模式解码（LoadIntoMemory=false）
%   OutDir         输出目录，默认 data/processed/level0
%   Plot           true=保存 ROI 预览 PNG，默认 true
%   ROIs           可选 {[r c h w], ...} 行/列/高/宽（相对读出的块坐标）
%
% 依赖：src/imaging/sentinel_l0_toolkit/（read_s1_raw.m 等）；
%       Python 解码环境在其原位 data/raw/Sentinel-L0条带/.venv。
% 注意：read_s1_raw 每次都会重新解码；整景解码耗时且写出 ~10GB MAT，请按需使用。
%       本脚本仅做“解码+复数读取+参数保存+ROI”，成像聚焦见 Read_Data_RD.m（第3周）。

% ---------- 解析参数 ----------
p = inputParser;
addParameter(p, 'RawRoot', '', @ischar);
addParameter(p, 'Polarization', 'VV', @(x)ismember(upper(x),{'HH','HV','VV','VH'}));
addParameter(p, 'StartPulse', 1, @(x)isscalar(x) && x>=1);
addParameter(p, 'NumPulses', 4000, @(x)isscalar(x) && x>=1);
addParameter(p, 'Full', false, @islogical);
addParameter(p, 'OutDir', 'D:/study/SAR/data/processed/level0', @ischar);
addParameter(p, 'Plot', true, @islogical);
addParameter(p, 'ROIs', {}, @iscell);
parse(p, varargin{:});
opt = p.Results;
pol = upper(opt.Polarization);

% ---------- 定位工具包 / Python 环境 ----------
toolkit = fullfile(fileparts(mfilename('fullpath')), 'sentinel_l0_toolkit');
addpath(toolkit);
assert(isfile(fullfile(toolkit,'read_s1_raw.m')), '未找到 %s/read_s1_raw.m', toolkit);
pyExe = locateVenv(toolkit);
fprintf('Python 解码环境：%s\n', pyExe);

% ---------- 定位主测量 .dat ----------
assert(~isempty(opt.RawRoot) && isfolder(opt.RawRoot), '请通过 RawRoot 指定 RAW 产品根目录');
meas = findMeasurementDat(opt.RawRoot, pol);
fprintf('主测量文件：%s\n', meas);

% ---------- 输出目录 ----------
mkdir(opt.OutDir);
mkdir(fullfile(opt.OutDir,'roi'));

% ---------- 解码 ----------
tag = sprintf('%s_start%d_n%d', lower(pol), opt.StartPulse, opt.NumPulses);
rawMat = fullfile(opt.OutDir, sprintf('level0_%s_raw.mat', tag));
if opt.Full
    [~, info] = read_s1_raw(meas, Inf, 'Polarization', pol, ...
        'Python', pyExe, 'OutputFile', rawMat, 'LoadIntoMemory', false, 'Plot', false);
    dims = size(matfile(rawMat), 'iq');
    params = buildParams(info, opt, meas, rawMat, '', dims, 'full');
    writeJSON(params, fullfile(opt.OutDir, sprintf('level0_%s_full_params.json', lower(pol))));
    fprintf('整景解码完成（磁盘模式）：%s\n', rawMat);
    return;
end

[echo, info] = read_s1_raw(meas, opt.NumPulses, 'Polarization', pol, ...
    'StartPulse', opt.StartPulse, 'Python', pyExe, ...
    'OutputFile', rawMat, 'LoadIntoMemory', true, 'Plot', false);
fprintf('已解码块：%d x %d（PRI x 距离样本）\n', size(echo,1), size(echo,2));

% 按 README 顺序：IQ 校正 -> SWST 距离窗对齐（仅内存矩阵可做）
[echo, info] = IQ_correction(echo, info);
[echo, info] = echo_SWST(echo, info);

% 保存校正后块
corrMat = fullfile(opt.OutDir, sprintf('level0_%s_corrected.mat', tag));
save(corrMat, 'echo', 'info', '-v7.3');

% 参数 JSON
dims = size(echo);
params = buildParams(info, opt, meas, rawMat, corrMat, dims, 'block');

% ---------- ROI ----------
rois = opt.ROIs;
if isempty(rois), rois = autoROIs(dims, 3); end
roiCell = {};
for k = 1:numel(rois)
    r = rois{k};
    r = [max(1,r(1)) max(1,r(2)) min(r(3),dims(1)-r(1)+1) min(r(4),dims(2)-r(2)+1)];
    sub = echo(r(1):r(1)+r(3)-1, r(2):r(2)+r(4)-1);
    save(fullfile(opt.OutDir,'roi',sprintf('roi_%d.mat',k)), 'sub', '-v7.3');
    if opt.Plot
        h = figure('Visible','off');
        imagesc(20*log10(max(abs(sub),1e-12)));
        axis xy; colormap gray; colorbar;
        title(sprintf('ROI %d (start_pulse=%d, col=%d)', k, r(1)+opt.StartPulse-1, r(2)));
        saveas(h, fullfile(opt.OutDir,'roi',sprintf('roi_%d.png',k)));
        close(h);
    end
    roiCell{k} = struct('start_pulse', r(1)+opt.StartPulse-1, ...
        'col', r(2), 'rows', r(3), 'cols', r(4));
end
params.rois = [roiCell{:}];

writeJSON(params, fullfile(opt.OutDir, sprintf('level0_%s_params.json', tag)));
fprintf('完成：%s\n  - 原始解码：%s\n  - 校正块  ：%s\n  - 参数     ：%s\n  - ROI     ：%s/roi\n', ...
    opt.OutDir, rawMat, corrMat, ...
    fullfile(opt.OutDir, sprintf('level0_%s_params.json', tag)), opt.OutDir);
end

function pyExe = locateVenv(toolkit)
cands = {
    fullfile(toolkit, '.venv', 'Scripts', 'python.exe'),
    'D:/study/SAR/data/raw/Sentinel-L0条带/.venv/Scripts/python.exe'
};
pyExe = '';
for i = 1:numel(cands)
    if isfile(cands{i}), pyExe = cands{i}; break; end
end
assert(~isempty(pyExe), '未找到 Python 解码环境（.venv）。请确认 data/raw/Sentinel-L0条带/.venv 存在。');
end

function meas = findMeasurementDat(rawRoot, pol)
pat = sprintf('s1?-s2-raw-s-%s-*.dat', lower(pol));
hits = dir(fullfile(rawRoot, '**', pat));
assert(numel(hits)>=1, '在 %s 下未找到 %s 极化的测量 .dat', rawRoot, pol);
keep = true(size(hits));
for i = 1:numel(hits)
    n = lower(hits(i).name);
    if contains(n,'-annot') || contains(n,'-index'), keep(i) = false; end
end
hits = hits(keep);
assert(numel(hits)>=1, '过滤 annot/index 后无主测量 .dat');
[~, idx] = max([hits.bytes]);
meas = fullfile(hits(idx).folder, hits(idx).name);
end

function params = buildParams(info, opt, meas, rawMat, corrMat, dims, mode)
r = info.radar;
params = struct();
params.product_id = 'S1A_S2_RAW__0SDV_20250303T140022_20250303T140047_058141_072E7E_7DD0';
params.sensor = 'Sentinel-1A';
params.mode = 'SM (stripmap) S2 beam';
params.polarization = upper(opt.Polarization);
params.measurement_file = meas;
params.decode_mode = mode;
params.carrier_frequency_hz = getf(r,'fc',[]);
params.wavelength_m = getf(r,'lambda',[]);
params.pri_s = getf(r,'pri_s',[]);
params.prf_hz = 1/getf(r,'pri_s',[]);
params.range_sample_rate_hz = getf(r,'range_sample_rate_hz',[]);
params.na_all = getf(info.image,'Na_all',[]);
params.rows = dims(1);
params.range_samples = dims(2);
params.start_pulse = opt.StartPulse;
params.num_pulses = opt.NumPulses;
params.output_mat_raw = rawMat;
if ~isempty(corrMat), params.output_mat_corrected = corrMat; end
params.generated_at = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
params.note = '第2周交付：Level-0 解码+复数读取+参数保存+3个ROI；整景聚焦见 Read_Data_RD.m（第3周）';
end

function v = getf(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function rois = autoROIs(sz, k)
nr = sz(1); nc = sz(2);
win = min(1024, nr-1); winc = min(1024, nc-1);
rois = {};
for i = 1:k
    t = (i-0.5)/k;
    r0 = max(1, floor((nr-win)*t)+1);
    c0 = max(1, floor((nc-winc)*t)+1);
    rois{end+1} = [r0, c0, win, winc];
end
end

function writeJSON(s, path)
json = jsonencode(s, 'PrettyPrint', true);
fid = fopen(path, 'w', 'n', 'UTF-8');
fwrite(fid, json);
fclose(fid);
end
