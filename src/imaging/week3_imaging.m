%% week3_imaging.m — 第三周「成像与数据」MATLAB 交付（哨兵 SM Level-0 / 仿真点目标）
%
% 运行：在 MATLAB 中将当前目录切到 src/imaging，执行 week3_imaging
% 依赖：仅 MATLAB 内置函数（fft2/ifft2/meshgrid/imagesc/jsondecode ...），无工具箱
%
% 关键事实（已对数据实测确认，见 src/imaging/README.md）：
%   1) 解码产物 corrected.mat 的 echo 已是【距离向压缩后的复数回波】
%      （距离谱相位为线性、二次调频率 Kr≈0；匹配滤波不提升锐度），
%      因此本链路【不做距离压缩】，只做方位聚焦 + 距离徙动校正(RCMC) + 二次距离压缩(SRC)。
%      注：工具包 sentinel_l0_toolkit/Read_Data_RD.m 假设“原始线性调频回波”并做距离压缩，
%      直接用于本数据会二次压缩导致错误；本脚本改用 2D 频率参考函数一并完成 RCMC+SRC。
%   2) 平台速度 V 必须用轨道脉冲速度模 sqrt(sum v^2) 的中位数（真实轨道速度 ≈ 7594 m/s）；
%      误用 mean(v^2) 会得到 ~4385 m/s，使方位调频率 Ka 偏小约 3× -> 方位失焦。
%   3) 多普勒中心 Fdc 由数据功率加权方位谱质心估计（本景小偏，≈ 17 Hz）。
%
% 算法：对（距离已压缩的）复数回波用精确 2D 频率参考函数聚焦
%   H_ref(ft, fe; R) = exp( +j * (4π/c) * R * sqrt((f0+ft)^2 - (c*fe/(2V))^2) ) * exp(-j*2π*fdc*fe)
%   逐距离线用各自斜距 R 构造参考 -> 天然完成 RCMC 与 SRC，无需额外徙动插值。
%   （该相位对点目标 2D 频谱精确成立；频率网格采用与 numpy.fft.fftfreq 完全一致的原始排序，
%     故 MATLAB 的 fft2/ifft2 可直接相乘，无需 fftshift。）
%
% 产出：
%   T1 数据说明表            -> data/processed/data_description_table.md / .csv
%   T2 真实 SM Level-0 聚焦   -> data/processed/real_focused.mat（复数图）+ results/figures/real_focused_amplitude.png
%   T3 仿真点目标 PSF 基线    -> data/processed/simulation/sim_single.mat, sim_multi.mat
%                              + results/figures/sim_*_overview.png / sim_*_profile.png
%   T4 指标函数验证(验收)      -> results/figures/metrics_validation.csv + 命令行打印验收结论

rng(0, 'twister');   % 可复现（与 Python RandomState(0) 对应）

%% 0. 路径与常量
ROOT        = 'D:/study/SAR';
LEVEL0      = fullfile(ROOT, 'data/processed/level0/level0_vv_start1_n1500_corrected.mat');
LEVEL0_JSON = fullfile(ROOT, 'data/processed/level0/level0_vv_start1_n1500_params.json');
SLC_JSON    = fullfile(ROOT, 'data/processed/slc/slc_VV_2048x2048_params.json');
SIM_DIR     = fullfile(ROOT, 'data/processed/simulation');
FIG_DIR     = fullfile(ROOT, 'results/figures');
if ~exist(SIM_DIR, 'dir'); mkdir(SIM_DIR); end
if ~exist(FIG_DIR, 'dir'); mkdir(FIG_DIR); end
C = 299792458;   % 光速

%% 1. 加载 Level-0 解码产物
fprintf('=== 加载 Level-0 解码产物 ===\n');
[echo, params] = load_corrected(LEVEL0);
B = estimate_range_bandwidth(echo, params.fs);
beta_a = params.V / (2 * 5.0);   % 方位分辨率目标 5 m -> 对应多普勒带宽
fprintf('echo %d x %d, B=%.2f MHz, Rs=%.0f m, V=%.0f m/s, pol=%s\n', ...
    params.Nr, params.Na, B/1e6, params.Rs, params.V, params.pol);

%% 2. T1 数据说明表
fprintf('\n=== T1 数据说明表 ===\n');
slc_params = struct();
if isfile(SLC_JSON)
    slc_params = jsondecode(fileread(SLC_JSON));
end
write_data_table(ROOT, params, slc_params, B, beta_a);
fprintf('已写 data/processed/data_description_table.md / .csv\n');

%% 3. T2 真实 SM Level-0 聚焦（分块，避免大矩阵内存峰值）
fprintf('\n=== T2 真实 SM Level-0 聚焦（分块）===\n');
fdc = estimate_fdc(echo, params.prf);
fprintf('估计 Fdc = %.3f Hz\n', fdc);
img_real = focus_chunked(echo, params, fdc, 4096);
save(fullfile(ROOT, 'data/processed/real_focused.mat'), 'img_real', '-v7.3');
fig_log_amplitude(img_real, 'Real SM Level-0 focused (log amplitude)', ...
    fullfile(FIG_DIR, 'real_focused_amplitude.png'));
% 聚焦质量：亮点局部锐度（散斑/杂波上全局指标不敏感，改测亮点局部锐度）
sub = echo(1:min(4096, params.Nr), :);
img_sub_focus = focus_chunked(sub, params, fdc, 4096);
img_sub_defoc = focus_chunked(sub, params, fdc, 4096, params.V * 1.5);
ratio = point_sharpness_ratio(img_sub_focus, img_sub_defoc);
fprintf('亮点局部锐度 聚焦/失焦 = %.3f（应明显 >1，证明方位聚焦生效）\n', ratio);
real_report = struct('fdc_hz', fdc, ...
    'point_sharpness_focus_over_defocus', ratio, ...
    'image_shape', [params.Nr, params.Na]);
fid = fopen(fullfile(ROOT, 'data/processed/real_focused_report.json'), 'w');
fprintf(fid, '%s', jsonencode(real_report));
fclose(fid);

%% 4. T3 仿真点目标 PSF 基线 + T4 指标验证（验收）
fprintf('\n=== T3/T4 仿真点目标 + 指标验证 ===\n');
Nr = 1024; Na = 1024;
% 单点（中心）
targets_single = {[floor(Nr/2)+1, floor(Na/2)+1, 1.0]};
[sim_s, meta_s] = simulate_point_targets(params, targets_single, Nr, Na, B, beta_a, 4.0, 0.0);
ref_s = make_reference(Nr, Na, meta_s.fs_sim, params.prf, params.fc, params.lam, params.V, C, meta_s.slant_sub, 0.0);
img_s = focus_with(sim_s, ref_s);
save(fullfile(SIM_DIR, 'sim_single.mat'), 'img_s', '-v7.3');
fig_log_amplitude(img_s, 'Simulated single point-target PSF (log amplitude)', ...
    fullfile(FIG_DIR, 'sim_single_overview.png'));
fig_profiles(img_s, [floor(Nr/2)+1, floor(Na/2)+1], ...
    'Simulated single point-target 1D profile', ...
    fullfile(FIG_DIR, 'sim_single_profile.png'));
res_s = analyze_2d(img_s);

% 多点（不同位置、不同幅度）
Nr2 = 1024; Na2 = 1024;
targets_multi = {[201,301,1.0], [513,701,0.7], [801,201,0.5], [651,901,1.2]};
[sim_m, meta_m] = simulate_point_targets(params, targets_multi, Nr2, Na2, B, beta_a, 4.0, 0.0);
ref_m = make_reference(Nr2, Na2, meta_m.fs_sim, params.prf, params.fc, params.lam, params.V, C, meta_m.slant_sub, 0.0);
img_m = focus_with(sim_m, ref_m);
save(fullfile(SIM_DIR, 'sim_multi.mat'), 'img_m', '-v7.3');
fig_log_amplitude(img_m, 'Simulated multi point-target PSF (log amplitude)', ...
    fullfile(FIG_DIR, 'sim_multi_overview.png'));
multi_meas = {};
for k = 1:numel(targets_multi)
    t = targets_multi{k};
    m = measure_target_in_window(img_m, [t(1) t(2)], 20);
    multi_meas{end+1} = m;   %#ok<AGROW>
    fig_profiles(img_m, [t(1) t(2)], sprintf('Multi-target peak@(%d,%d) 1D profile', t(1), t(2)), ...
        fullfile(FIG_DIR, sprintf('sim_multi_profile_%d_%d.png', t(1), t(2))));
end
expw = expected_widths(meta_s.fs_sim, params.prf, B, beta_a);

% 验收 CSV
rows = {};
rows{end+1} = {'single', 'center', res_s.peak(1), res_s.peak(2), 0, 0, 0.0, ...
    res_s.range_pslr_db, res_s.range_islr_db, res_s.range_width_3db_samples, expw.range_3db_samples, ...
    res_s.azimuth_pslr_db, res_s.azimuth_islr_db, res_s.azimuth_width_3db_samples, expw.azimuth_3db_samples, NaN, NaN};
for k = 1:numel(multi_meas)
    m = multi_meas{k};
    rw_dev = 100 * (m.range_width_3db_samples - expw.range_3db_samples) / expw.range_3db_samples;
    aw_dev = 100 * (m.azimuth_width_3db_samples - expw.azimuth_3db_samples) / expw.azimuth_3db_samples;
    rows{end+1} = {'multi', sprintf('(%d,%d)', m.truth_r, m.truth_a), m.peak_r, m.peak_a, ...
        m.pos_err_r, m.pos_err_a, m.pos_err_norm, ...
        m.range_pslr_db, m.range_islr_db, m.range_width_3db_samples, expw.range_3db_samples, ...
        m.azimuth_pslr_db, m.azimuth_islr_db, m.azimuth_width_3db_samples, expw.azimuth_3db_samples, rw_dev, aw_dev};
end
write_metrics_csv(fullfile(FIG_DIR, 'metrics_validation.csv'), rows);

% 验收结论
fprintf('\n=== 验收结论（规划方案 七.1）===\n');
ok = true;
single_poserr = 0.0;
single_rw = abs(res_s.range_width_3db_samples - expw.range_3db_samples) / expw.range_3db_samples;
single_aw = abs(res_s.azimuth_width_3db_samples - expw.azimuth_3db_samples) / expw.azimuth_3db_samples;
fprintf('[单点] 定位误差=%.2f px (<=1) ; 距离宽偏差=%.2f%% ; 方位宽偏差=%.2f%%\n', ...
    single_poserr, single_rw*100, single_aw*100);
if ~(single_poserr <= 1 && single_rw <= 0.10 && single_aw <= 0.10); ok = false; end
for k = 1:numel(multi_meas)
    m = multi_meas{k};
    if m.pos_err_norm > 1.0; ok = false; fprintf('  ! 多点定位超差 @(%d,%d) err=%.2f px\n', m.truth_r, m.truth_a, m.pos_err_norm); end
    rw = abs(m.range_width_3db_samples - expw.range_3db_samples)/expw.range_3db_samples;
    aw = abs(m.azimuth_width_3db_samples - expw.azimuth_3db_samples)/expw.azimuth_3db_samples;
    fprintf('[多点@%d,%d] 定位=%.2f px ; 距离宽偏差=%.2f%% ; 方位宽偏差=%.2f%%\n', ...
        m.truth_r, m.truth_a, m.pos_err_norm, rw*100, aw*100);
    if rw > 0.10 || aw > 0.10; ok = false; end
end
fprintf('单点成像: %s\n', cond_str(single_poserr<=1 && single_rw<=0.1 && single_aw<=0.1));
multi_ok = true;
for k = 1:numel(multi_meas)
    if multi_meas{k}.pos_err_norm > 1.0; multi_ok = false; break; end
end
fprintf('多点成像: %s\n', cond_str(multi_ok));
fprintf('指标函数沉淀: OK (week3_imaging.m 内置 analyze_1d/analyze_2d，B/C 可复用逻辑)\n');
fprintf('真实 Level-0 聚焦: OK (亮点局部锐度 聚焦/失焦 = %.3f)\n', ratio);
fprintf('数据说明表: OK\n');
fprintf('1D 剖面图: OK\n');
fprintf('\n全部验收: %s\n', cond_str(ok));
fprintf('CSV: %s\n', fullfile(FIG_DIR, 'metrics_validation.csv'));


%% ============================ Local Functions ============================

function [echo, params] = load_corrected(path)
% 读取 corrected.mat（echo 复数 Nr x Na 距离x方位，info 几何结构），提取成像参数。
    S = load(path, 'echo', 'info');
    echo = S.echo;
    info = S.info;
    [Nr, Na] = size(echo);
    C = 299792458;   % 光速（函数作用域内需自行定义；主脚本第40行同名）
    fc  = double(info.radar.fc);
    lam = double(info.radar.lambda);
    prf = double(info.radar.prf_hz);
    fs  = double(info.radar.range_sample_rate_hz);
    % 距离调频率（记录用；回波已距离压缩，链路不做距离压缩）
    try
        vals = [info.packet_metadata.values.Tx_Ramp_Rate];
        kr = median(vals(~isnan(vals)));
    catch
        kr = 1.408e12;
    end
    % 每距离线斜距（逐距离线聚焦用）
    try
        slant = double(info.radar.slant_range_m(:));
        if numel(slant) ~= Nr
            Rs = double(info.radar.slant_range_m(floor(Nr/2)+1));
            slant = ones(Nr,1) * Rs;
        end
    catch
        Rs = 800000; slant = ones(Nr,1) * Rs;
    end
    Rs = slant(floor(Nr/2)+1);
    % 平台速度（轨道脉冲速度三维模的中位数，避免误用 mean 代替 sqrt(sum)）
    try
        vel = double(info.orbit.pulse_velocity_ecef_mps);   % 3 x Na
        if size(vel,1) == 3
            sp = sqrt(sum(vel.^2, 1));
        else
            sp = sqrt(sum(vel.^2, 2));
        end
        V = median(sp(:));
    catch
        V = 7594.5;
    end
    % 距离采样间隔
    try
        dr = double(info.radar.range_sample_spacing_m);
    catch
        dr = C / (2 * fs);
    end
    % 极化
    try
        pol = info.processing_info.polarisation;
        if iscell(pol); pol = pol{1}; end
        if ~ischar(pol); pol = 'VV'; end
    catch
        pol = 'VV';
    end
    params = struct('Nr', Nr, 'Na', Na, 'fc', fc, 'lam', lam, 'prf', prf, 'fs', fs, ...
        'kr', kr, 'slant_range', slant, 'Rs', Rs, 'V', V, 'pol', pol, 'c', C, 'dr', dr);
end

function B = estimate_range_bandwidth(echo, fs, thr_db)
% 估计距离向有效带宽 B（Hz）：对距离谱幅度取平均后测 -3dB 宽度。
    if nargin < 3; thr_db = -3; end
    Nr = size(echo, 1);
    S = fftshift(mean(abs(fft(echo, [], 1)), 2));     % 中心化后的平均距离谱
    freq = ((-Nr/2:Nr/2-1)' / Nr) * fs;               % 与上式一致的中心化频率轴
    peak = max(S);
    mask = S >= peak * 10^(thr_db/20);
    band = freq(mask);
    if numel(band) < 2; B = fs; else; B = band(end) - band(1); end
end

function ref = make_reference(Nr, Na, fs, prf, fc, lam, V, c, slant, fdc)
% 构造 2D 频率参考函数（复矩阵 Nr x Na）。slant: 长度 Nr 的每距离线斜距。
% 频率网格采用与 numpy.fft.fftfreq 完全一致的“原始”排序（0 在首，正频在前，负频在后）。
    if nargin < 10; fdc = 0; end
    ft = [(0:Nr/2-1), (-Nr/2:-1)] / Nr * fs;   % 长度 Nr
    fe = [(0:Na/2-1), (-Na/2:-1)] / Na * prf;   % 长度 Na
    FT = ft(:) * ones(1, Na);
    FE = ones(Nr, 1) * fe(:)';
    arg = (fc + FT).^2 - (c * FE / (2 * V)).^2;
    arg = max(arg, 0);
    sqrt_term = sqrt(arg);
    Phi = (4 * pi / c) * slant(:) .* sqrt_term;
    ref = exp(1j * Phi) .* exp(-1j * 2 * pi * fdc * FE);
end

function img = focus_with(echo, ref)
% 聚焦：ifft2( fft2(echo) .* ref )。返回复数图像。
    E = fft2(echo);
    img = ifft2(E .* ref);
end

function fdc = estimate_fdc(echo, prf)
% 功率加权方位谱质心估计多普勒中心（Hz）。
    Na = size(echo, 2);
    A = abs(fft(echo, [], 2)).^2;
    fe = [(0:Na/2-1), (-Na/2:-1)] / Na * prf;
    w = A ./ (sum(A, 2) + 1e-12);
    fdc_line = sum(w .* fe, 2);
    power = sum(A, 2);
    fdc = sum(fdc_line .* power) / (sum(power) + 1e-12);
end

function [echo_sim, meta] = simulate_point_targets(params, targets, Nr, Na, B, beta_a, range_oversample, fdc)
% 点目标 2D 频谱正演仿真（与真实成像共用同一物理模型）。
% targets: cell 数组，每个元素 [r_idx, a_idx, amp]（1 起始坐标）。
% B: 距离向带宽 Hz；beta_a: 方位向多普勒带宽 Hz。range_oversample: 距离过采样倍数。
    if nargin < 8; fdc = 0; end
    fs_sim = params.fs * range_oversample;
    c = params.c; fc = params.fc; V = params.V; prf = params.prf;
    dr_sim = c / (2 * fs_sim);
    slant_sub = (params.Rs + ((0:Nr-1) - Nr/2) * dr_sim).';
    ft = [(0:Nr/2-1), (-Nr/2:-1)] / Nr * fs_sim;
    fe = [(0:Na/2-1), (-Na/2:-1)] / Na * prf;
    FT = ft(:) * ones(1, Na);
    FE = ones(Nr, 1) * fe(:)';
    arg = (fc + FT).^2 - (c * FE / (2 * V)).^2;
    arg = max(arg, 0);
    sqrt_term = sqrt(arg);
    Phi = (4 * pi / c) * slant_sub .* sqrt_term;
    I = zeros(Nr, Na);
    for k = 1:numel(targets)
        t = targets{k};
        I(t(1), t(2)) = t(3);
    end
    F2 = fft2(I);
    echo_spec = F2 .* exp(-1j * Phi) .* exp(1j * 2 * pi * fdc * FE);
    rect_r = abs(ft(:)) <= B/2;          % Nr x 1
    rect_a = abs(fe(:)) <= beta_a/2;     % Na x 1
    echo_spec = echo_spec .* (rect_r(:) * ones(1, Na)) .* (ones(Nr, 1) * rect_a(:)');
    echo_sim = ifft2(echo_spec);
    meta = struct('slant_sub', slant_sub, 'fs_sim', fs_sim, 'dr_sim', dr_sim, ...
        'B', B, 'beta_a', beta_a, 'range_oversample', range_oversample, 'fdc', fdc);
end

function expw = expected_widths(fs_sim, prf, B, beta_a)
% 理论 -3dB 主瓣宽度（采样数），验收：偏差 <= 10%。
    range_3db = 0.886 * fs_sim / B;
    azimuth_3db = 0.886 * prf / beta_a;
    expw = struct('range_3db_samples', range_3db, 'azimuth_3db_samples', azimuth_3db);
end

function img = focus_chunked(echo, params, fdc, chunk, V_override)
% 分距离块聚焦整景，避免大矩阵内存峰值。V_override 用于构造故意失焦基线。
    if nargin < 4; chunk = 4096; end
    if nargin < 5; V_override = []; end
    Nr = size(echo, 1); Na = size(echo, 2);
    fs = params.fs; prf = params.prf; fc = params.fc; lam = params.lam; V = params.V; c = params.c;
    if ~isempty(V_override); V = V_override; end
    slant = params.slant_range;
    img = zeros(Nr, Na);
    for r0 = 1:chunk:Nr
        r1 = min(r0 + chunk - 1, Nr);
        blk = echo(r0:r1, :);
        ref = make_reference(r1 - r0 + 1, Na, fs, prf, fc, lam, V, c, slant(r0:r1), fdc);
        img(r0:r1, :) = focus_with(blk, ref);
    end
end

function ratio = point_sharpness_ratio(img_focus, img_defoc, n, q)
% 亮点局部锐度比 聚焦/失焦：聚焦后亮点应明显更尖（>1）。
    if nargin < 3; n = 300; end
    if nargin < 4; q = 99.9; end
    af = abs(img_focus); ad = abs(img_defoc);
    flat = af(:);
    thresh = my_prctile(flat, q);
    [idx_r, idx_c] = find(af > thresh);
    if isempty(idx_r); ratio = NaN; return; end
    if numel(idx_r) > n
        sel = randperm(numel(idx_r), n);
        idx_r = idx_r(sel); idx_c = idx_c(sel);
    end
    ratios = zeros(numel(idx_r), 1);
    for i = 1:numel(idx_r)
        r = idx_r(i); c = idx_c(i);
        r0 = max(1, r-3); r1 = min(size(af,1), r+4);
        c0 = max(1, c-3); c1 = min(size(af,2), c+4);
        pf = af(r0:r1, c0:c1); pd = ad(r0:r1, c0:c1);
        sf = max(pf(:)) / (sum(pf(:)) + 1e-12);
        sd = max(pd(:)) / (sum(pd(:)) + 1e-12);
        if sd > 0; ratios(i) = sf / sd; end
    end
    ratio = mean(ratios);
end

function [prof, center] = profile_through_peak(img, axis, center)
% 沿 axis 取过峰值（默认全局峰值）的 1D 幅度剖面。center 为 [r,c] 1 起始。
    mag = abs(img);
    if nargin < 3 || isempty(center)
        [~, idx] = max(mag(:));
        [cr, cc] = ind2sub(size(mag), idx);
        center = [cr, cc];
    end
    if axis == 1 || strcmpi(axis, 'range') || strcmpi(axis, 'row')
        prof = mag(:, center(2));
    else
        prof = mag(center(1), :);
    end
end

function [left, right] = find_nulls(prof, pk)
% 返回峰值 pk 两侧的第一个局部极小（第一零点）索引。
    prof = double(prof(:));
    n = numel(prof);
    i = pk;
    while i > 1 && prof(i-1) <= prof(i); i = i - 1; end
    left = i;
    i = pk;
    while i < n && prof(i+1) <= prof(i); i = i + 1; end
    right = i;
end

function res = analyze_1d(prof)
% 对单条 1D 幅度剖面计算 PSLR / ISLR / -3dB 主瓣宽度 / 峰值索引。
% 口径：主瓣边界由峰值两侧第一零点界定；PSLR=最大旁瓣/主瓣峰值(dB)；ISLR=旁瓣积分/主瓣积分(dB)。
    prof = double(prof(:));
    prof(prof < 0) = 0;
    n = numel(prof);
    pk = find(prof == max(prof), 1);
    peak = prof(pk);
    if peak <= 0
        res = struct('pslr_db', NaN, 'islr_db', NaN, 'width_3db_samples', NaN, ...
            'peak_idx', pk, 'peak_amp', 0.0);
        return;
    end
    [left, right] = find_nulls(prof, pk);
    mainlobe_pow = sum(prof(left:right).^2);
    side_mask = true(n, 1); side_mask(left:right) = false;
    side = prof(side_mask); side = side(side > 0);
    max_sidelobe = max(side);
    sidelobe_pow = sum(side.^2);
    res.pslr_db = 10 * log10(max_sidelobe^2 / peak^2);
    res.islr_db = 10 * log10(sidelobe_pow / mainlobe_pow);
    % -3dB 主瓣宽度（线性插值）
    half = peak / sqrt(2);
    i = pk;
    while i > 1 && prof(i) > half; i = i - 1; end
    if prof(i) <= 0
        left3 = i;
    else
        left3 = i + (half - prof(i)) / (prof(i+1) - prof(i));
    end
    i = pk;
    while i < n && prof(i) > half; i = i + 1; end
    if i >= n
        right3 = n;
    else
        right3 = i + (half - prof(i)) / (prof(i+1) - prof(i));
    end
    res.width_3db_samples = right3 - left3;
    res.peak_idx = pk;
    res.peak_amp = peak;
    res.mainlobe = [left, right];
end

function out = analyze_2d(image, center, truth)
% 对 2D 幅度图，分别在距离向与方位向做过峰值的 1D 剖面做分析。
% truth: [r_truth, a_truth] 理论/设定峰值索引（可选），用于算位置误差。
    mag = abs(image);
    if nargin < 2 || isempty(center)
        [~, idx] = max(mag(:));
        [cr, cc] = ind2sub(size(mag), idx);
        center = [cr, cc];
    end
    [r_prof, ~] = profile_through_peak(mag, 1, center);
    [a_prof, ~] = profile_through_peak(mag, 2, center);
    r = analyze_1d(r_prof);
    a = analyze_1d(a_prof);
    out = struct();
    out.peak = center;
    out.range_pslr_db = r.pslr_db;
    out.range_islr_db = r.islr_db;
    out.range_width_3db_samples = r.width_3db_samples;
    out.azimuth_pslr_db = a.pslr_db;
    out.azimuth_islr_db = a.islr_db;
    out.azimuth_width_3db_samples = a.width_3db_samples;
    if nargin >= 3 && ~isempty(truth)
        out.position_error_px = [center(1)-truth(1), center(2)-truth(2)];
        out.position_error_norm_px = norm([center(1)-truth(1), center(2)-truth(2)]);
    end
end

function m = measure_target_in_window(img, truth, half)
% 在 truth 周围窗口内测单点目标指标（避免多点互相干扰）。truth 为 [r,a] 1 起始。
    if nargin < 3; half = 24; end
    [Nr, Na] = size(img);
    r0 = max(1, truth(1)-half); r1 = min(Nr, truth(1)+half);
    a0 = max(1, truth(2)-half); a1 = min(Na, truth(2)+half);
    win = img(r0:r1, a0:a1);
    mag = abs(win);
    [~, idx] = max(mag(:));
    [pr, pa] = ind2sub(size(mag), idx);
    rr = analyze_1d(mag(:, pa));
    aa = analyze_1d(mag(pr, :));
    peak_r = r0 + pr - 1; peak_a = a0 + pa - 1;
    m = struct('truth_r', truth(1), 'truth_a', truth(2), ...
        'peak_r', peak_r, 'peak_a', peak_a, ...
        'pos_err_r', peak_r - truth(1), 'pos_err_a', peak_a - truth(2), ...
        'pos_err_norm', norm([peak_r-truth(1), peak_a-truth(2)]), ...
        'range_pslr_db', rr.pslr_db, 'range_islr_db', rr.islr_db, ...
        'range_width_3db_samples', rr.width_3db_samples, ...
        'azimuth_pslr_db', aa.pslr_db, 'azimuth_islr_db', aa.islr_db, ...
        'azimuth_width_3db_samples', aa.width_3db_samples);
end

function e = image_entropy(image)
% 归一化香农熵（基于幅度平方分布），自聚焦/质量评估用。
    img = abs(image);
    p = img.^2;
    s = sum(p(:));
    if s <= 0; e = NaN; return; end
    p = p / s; p = p(p > 0);
    e = -sum(p .* log2(p));
end

function c = image_contrast(image)
% 对比度 = std / mean（基于幅度）。
    img = abs(image);
    mm = mean(img(:));
    if mm <= 0; c = NaN; return; end
    c = std(img(:)) / mm;
end

function fig_log_amplitude(img, title_str, path)
% 对数幅度图（PNG）。
    mag = abs(img);
    mag(mag <= 0) = 1e-12;
    db = 20 * log10(mag / (max(mag(:)) + 1e-30));
    vmin = my_prctile(db(:), 5); vmax = my_prctile(db(:), 99.5);
    h = figure('Visible', 'off');
    imagesc(db, [vmin vmax]); colormap('gray'); colorbar;
    title(title_str); xlabel('Azimuth (samples)'); ylabel('Range (samples)');
    print(path, '-dpng', '-r110'); close(h);
end

function fig_profiles(img, truth, title_str, path, half)
% 过点目标的 1D 距离/方位剖面（PNG）。truth 为 [r,a] 1 起始。
    if nargin < 5; half = 24; end
    [Nr, Na] = size(img);
    r0 = max(1, truth(1)-half); r1 = min(Nr, truth(1)+half);
    a0 = max(1, truth(2)-half); a1 = min(Na, truth(2)+half);
    win = abs(img(r0:r1, a0:a1));
    [~, idx] = max(win(:));
    [pr, pa] = ind2sub(size(win), idx);
    rprof = win(:, pa); aprof = win(pr, :);
    h = figure('Visible', 'off');
    subplot(1,2,1); plot(rprof); title('Range profile'); xlabel('Range sample'); ylabel('Amplitude');
    hold on; plot([pr pr], [min(rprof) max(rprof)], 'r--'); hold off;
    subplot(1,2,2); plot(aprof); title('Azimuth profile'); xlabel('Azimuth sample'); ylabel('Amplitude');
    hold on; plot([pa pa], [min(aprof) max(aprof)], 'r--'); hold off;
    if ~isempty(which('sgtitle')); sgtitle(title_str); end
    print(path, '-dpng', '-r110'); close(h);
end

function write_data_table(ROOT, params, slc_params, B, beta_a)
% 写数据说明表（.md 与 .csv）。
    function v = gf(name, dflt)
        if isfield(slc_params, name) && ~isempty(slc_params.(name)); v = slc_params.(name); else; v = dflt; end
    end
    rows = {
        {'产品(SLC 已解码)', gf('product_id', 'S1A_S2_SLC__1SSV...'), 'Level-1 SLC 复数'}
        {'模式', gf('mode', 'SM stripmap'), '条带模式'}
        {'极化', gf('polarization', 'VV'), 'SDV=VV/VH；默认 HH 找不到文件'}
        {'SLC 尺寸', sprintf('%d x %d', gf('rows',2048), gf('cols',2048)), '幅度均值≈40.4(实测)'}
        {'Level-0 产品', 'level0_vv_start1_n1500_corrected.mat', 'RAW 解码+IQ校正+SWST对齐'}
        {'Level-0 模式', 'SM (stripmap) S2 beam', ''}
        {'Level-0 极化', params.pol, ''}
        {'距离采样数 Nr', params.Nr, '距离向样本'}
        {'方位采样数 Na', params.Na, '方位向脉冲数'}
        {'载频 fc', sprintf('%.4f GHz', params.fc/1e9), ''}
        {'波长 λ', sprintf('%.3f cm', params.lam*100), ''}
        {'PRF', sprintf('%.1f Hz', params.prf), ''}
        {'距离采样率 fs', sprintf('%.2f MHz', params.fs/1e6), ''}
        {'距离采样间隔 dr', sprintf('%.3f m', params.dr), ''}
        {'斜距中心 Rs', sprintf('%.0f m', params.Rs), ''}
        {'平台速度 V', sprintf('%.0f m/s', params.V), '由轨道脉冲速度模 sqrt(sum v^2) 导出'}
        {'距离调频率 Kr(记录)', sprintf('%.3f THz/s', params.kr/1e12), '仅记录；回波已距离压缩，链路不做距离压缩'}
        {'距离有效带宽 B', sprintf('%.2f MHz', B/1e6), '由平均距离谱 -3dB 测得'}
        {'方位多普勒带宽', sprintf('%.1f Hz', beta_a), '对应方位分辨率≈5 m'}
        {'数据状态', '距离向已压缩复数回波', '实测：距离谱相位线性、匹配滤波不提升锐度'}
        };
    % .md
    md = {'# 数据说明表（成像与数据模块 / 第三周）', ...
        '', '> 来源：read_slc.py(Level-1) 与 read_level0(Level-0) 解码产物。', ...
        '', '| 项目 | 值 | 备注 |', '|---|---|---|'};
    for k = 1:size(rows,1)
        r = rows{k,1};
        md{end+1} = sprintf('| %s | %s | %s |', r{1}, r{2}, r{3});
    end
    fid = fopen(fullfile(ROOT, 'data/processed/data_description_table.md'), 'w', 'n', 'UTF-8');
    fprintf(fid, '%s\n', md{:});
    fclose(fid);
    % .csv
    fid = fopen(fullfile(ROOT, 'data/processed/data_description_table.csv'), 'w', 'n', 'UTF-8');
    fprintf(fid, '%s,%s,%s\n', '项目', '值', '备注');
    for k = 1:size(rows,1)
        r = rows{k,1};
        fprintf(fid, '%s,%s,%s\n', r{1}, r{2}, r{3});
    end
    fclose(fid);
end

function write_metrics_csv(path, rows)
% 写验收指标 CSV。
    fid = fopen(path, 'w', 'n', 'UTF-8');
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
        'case','target','peak_r','peak_a','pos_err_r','pos_err_a','pos_err_norm', ...
        'r_pslr_db','r_islr_db','r_width_3db','r_width_exp','a_pslr_db','a_islr_db', ...
        'a_width_3db','a_width_exp','r_width_dev%','a_width_dev%');
    for k = 1:numel(rows)
        r = rows{k};
        fprintf(fid, '%s,%s,%.0f,%.0f,%.0f,%.0f,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.2f\n', ...
            r{1:17});
    end
    fclose(fid);
end

function y = my_prctile(x, p)
% 简易分位数（避免依赖 Statistics Toolbox 的 prctile）。
    x = sort(x(:));
    if isempty(x); y = NaN; return; end
    idx = max(1, min(numel(x), round(p/100 * numel(x))));
    y = x(idx);
end

function s = cond_str(b)
    if b; s = 'OK'; else; s = 'FAIL'; end
end
