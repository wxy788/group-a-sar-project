%% week3_autofocus.m — 第三周「自聚焦」MATLAB 交付（散焦注入 + 最小熵自聚焦 + PGA 基线）
%
% 运行：在 MATLAB 中将当前目录切到 src/autofocus，执行 week3_autofocus
% 依赖：仅 MATLAB 内置函数（fft/ifft/imagesc/fminsearch/optimset ...），无工具箱
%
% 对应规划方案（哨兵SAR图像处理项目规划方案·缩减版）：
%   §二.2  阶段二 深度学习自聚焦 —— 第 3 周要求「散焦注入 + PGA 完成，建立剖面和指标」
%   §三.3  样本构造：使用条带仿真聚焦图和真实 SM SLC 复数小块，在方位频谱中注入已知
%          相位误差，保存散焦输入、参考图及误差参数
%   §三.4  指标：幅度 NRMSE、图像熵、图像对比度、点目标主瓣宽度、耗时
%   §七.1  验收：至少 6 个独立非零误差案例，其中至少 5 个 NRMSE 优于散焦输入；
%          另测零误差和未见误差；与传统方法比较质量及耗时
%
% 物理模型（与 PAFnet 等文献一致，误差加在「方位频谱」上）：
%   散焦  D = ifft( fft(I0, [], 2) .* exp(+j*phi(f)) , [], 2)
%   补偿  I = ifft( fft(D , [], 2) .* exp(-j*phi_hat(f)) , [], 2)
%   其中 f 为方位频率（列方向即方位向），phi 为待注入的已知相位误差。
%   该模型等价于「真实图像与散焦核 h=ifft(exp(j*phi)) 卷积」，是自聚焦的标准可逆模型：
%   只要 phi 已知即可精确复原，因此可用幅度 NRMSE 作为「恢复能力」的客观度量。
%
% 两条传统基线：
%   1) 最小熵自聚焦（MEAF）：phi_hat = B*c，B 为「偶次多项式 + 余弦(DCT)」基，
%      以图像熵为目标做 Nelder-Mead 无约束极小化（先粗扫二次项再进 fminsearch）。
%      优点：不需要场景中有强孤立点目标；缺点：受基函数表达能力限制。
%   2) PGA（相位梯度自聚焦）：选能量最强的若干距离单元，把每行亮点圆移到 index 1 对齐，
%      加环绕升余弦窗，在方位谱上取相位梯度 Δφ(k)=angle(Σ X(k)·conj(X(k-1)))，带内累积。
%      优点：非参数、可捕捉高频误差；缺点：依赖强孤立散射体，对分布杂波会失效（故加
%      「熵不劣化守卫」防止把已聚焦图改坏）。
%
%      ⚠ 两个实测踩过的坑（已在实现中修正，勿回退）：
%      · 峰值必须圆移到 index 1（不是数组中心）。移到中心会引入 (-1)^k 的 π 交替相位，
%        相位梯度出现系统性 ±π 跳变，累积后图像被彻底破坏（NRMSE≈1.4）。
%      · unwrap 与去线性必须【只在有效频带内】做。方位谱是带限的（仿真/实测均约 47%），
%        带外 |X|≈0、相位纯噪声；若在全带 unwrap/去趋势，会得到错尺度甚至反向的相位
%        （实测估计值只有真值的 ~0.5 倍且整体崩坏）。
%
% 产出：
%   data/processed/autofocus/cases/af_case_*.mat        9 组案例的 参考/散焦/复原/phi
%   data/processed/autofocus/samples/*.mat              第4周网络训练用 256×256 散焦-参考对
%   data/processed/autofocus/dataset_index.csv          样本索引（误差类型/参数/来源）
%   results/metrics/autofocus_metrics.csv               逐案例指标表
%   results/figures/autofocus_summary.png               熵/NRMSE 汇总
%   results/figures/autofocus_case_*.png                代表性案例 2×2 对比图
%   results/figures/autofocus_entropy_curves.png        PGA 迭代熵曲线
%   data/processed/autofocus/autofocus_week3_report.json 验收结论

rng(0, 'twister');   % 可复现

%% 0. 路径与常量
ROOT       = 'D:/study/SAR';
SLC_MAT    = fullfile(ROOT, 'data/processed/slc/slc_VV_2048x2048.mat');
SIM_MAT    = fullfile(ROOT, 'data/processed/simulation/sim_multi.mat');
FOCUS_MAT  = fullfile(ROOT, 'data/processed/real_focused.mat');
AF_DIR     = fullfile(ROOT, 'data/processed/autofocus');
CASE_DIR   = fullfile(AF_DIR, 'cases');
SAMP_DIR   = fullfile(AF_DIR, 'samples');
FIG_DIR    = fullfile(ROOT, 'results/figures');
MET_DIR    = fullfile(ROOT, 'results/metrics');
for d = {AF_DIR, CASE_DIR, SAMP_DIR, FIG_DIR, MET_DIR}
    if ~exist(d{1}, 'dir'); mkdir(d{1}); end
end

BLK = 512;              % 自聚焦处理块尺寸（距离 × 方位）
Na  = BLK;

cfg_meaf = struct('maxfev', 2500, 'coarse_n', 33, 'coarse_max', 16);
cfg_pga  = struct('bin_frac', 0.15, 'min_bins', 24, 'win_frac', [0.25 0.1667 0.125 0.125 0.10], ...
                  'band_thr', 0.05);

fprintf('=== 第三周 自聚焦基线（散焦注入 + 最小熵 + PGA）===\n');
fprintf('块尺寸 %d x %d（距离 x 方位），MEAF 最大评估次数 %d，PGA 迭代 %d 次\n\n', ...
    BLK, BLK, cfg_meaf.maxfev, numel(cfg_pga.win_frac));

%% 1. 三张底图（参考聚焦像）
fprintf('=== 1. 构造底图（参考聚焦像）===\n');
bases = build_bases(SIM_MAT, SLC_MAT, FOCUS_MAT, BLK);
names = fieldnames(bases);
for k = 1:numel(names)
    b = bases.(names{k});
    fprintf('  %-6s %d x %d  来源=%s  尖峰度(max/median)=%.2f\n', ...
        names{k}, size(b.img,1), size(b.img,2), b.src, b.spikiness);
end

%% 2. 误差案例定义（8 个非零 + 1 个零误差 + 1 个未见误差）
fprintf('\n=== 2. 相位误差案例 ===\n');
cases = define_cases();
for k = 1:numel(cases)
    fprintf('  %-3s base=%-3s type=%-9s %s\n', cases(k).id, cases(k).base, cases(k).type, cases(k).desc);
end

%% 3. 逐案例：注入 -> MEAF -> PGA -> 指标
fprintf('\n=== 3. 散焦注入与自聚焦 ===\n');
nC = numel(cases);
R  = struct([]);
for k = 1:nC
    cs   = cases(k);
    ref  = bases.(cs.base).img;
    u    = azimuth_u(Na);
    phi  = make_error(u, cs, cs.seed);
    def  = apply_phase(ref, phi);         % 方位频谱乘 exp(+j*phi) 注入已知误差

    fprintf('[%d/%d] %-3s %-9s 注入幅度 RMS=%.3f rad 峰值=%.3f rad ... ', ...
        k, nC, cs.id, cs.type, sqrt(mean(phi.^2)), max(abs(phi)));

    t0 = tic;
    [img_meaf, phi_meaf, ~] = meaf_autofocus(def, u, cfg_meaf);
    t_meaf = toc(t0);

    t0 = tic;
    [img_pga, phi_pga, ehist_pga, n_acc] = pga_autofocus(def, cfg_pga);
    t_pga = toc(t0);
    if n_acc < 1
        % PGA 守卫判定「补偿会更差」——保持散焦输入（熵指标下无改善空间）
        ehist_pga = image_entropy(def);
    else
        ehist_pga = [image_entropy(def), ehist_pga(1:n_acc)];
    end

    r = struct();
    r.id = cs.id;  r.base = cs.base;  r.type = cs.type;  r.desc = cs.desc;
    r.phi_rms_inj = sqrt(mean(phi.^2));
    r.entropy_ref    = image_entropy(ref);
    r.entropy_def    = image_entropy(def);
    r.entropy_meaf   = image_entropy(img_meaf);
    r.entropy_pga    = image_entropy(img_pga);
    r.contrast_ref   = image_contrast(ref);
    r.contrast_def   = image_contrast(def);
    r.contrast_meaf  = image_contrast(img_meaf);
    r.contrast_pga   = image_contrast(img_pga);
    r.nrmse_def      = nrmse_amp(ref, def);
    r.nrmse_meaf     = nrmse_amp(ref, img_meaf);
    r.nrmse_pga      = nrmse_amp(ref, img_pga);
    r.time_meaf_s    = t_meaf;
    r.time_pga_s     = t_pga;
    r.phi_resid_meaf = phase_resid_rms(phi_meaf, phi);
    r.phi_resid_pga  = phase_resid_rms(phi_pga,  phi);
    [r.width_ref, r.width_def, r.width_meaf, r.width_pga] = ...
        azimuth_width3db(ref, def, img_meaf, img_pga);
    r.entropy_hist_pga = ehist_pga;
    if isempty(R); R = r; else; R(end+1) = r; end   %#ok<SAGROW>

    fprintf('NRMSE 散焦=%.4f -> MEAF=%.4f (%.2fs) / PGA=%.4f (%.2fs)\n', ...
        r.nrmse_def, r.nrmse_meaf, t_meaf, r.nrmse_pga, t_pga);

    save(fullfile(CASE_DIR, sprintf('af_case_%s.mat', cs.id)), ...
        'ref', 'def', 'img_meaf', 'img_pga', 'phi', 'phi_meaf', 'phi_pga', ...
        'cs', 'r', '-v7.3');
end

%% 4. 第 4 周网络训练用散焦样本（256×256，真实 SLC）
fprintf('\n=== 4. 生成网络训练用散焦-参考样本 ===\n');
samp_info = make_training_samples(SLC_MAT, SAMP_DIR, 256, 12, 2);
fprintf('  已生成 %d 对 256x256 复数样本 -> %s\n', numel(samp_info), SAMP_DIR);

%% 5. 写出指标 CSV / 索引 CSV / JSON / 图
fprintf('\n=== 5. 输出指标与图 ===\n');
write_metrics_csv(fullfile(MET_DIR, 'autofocus_metrics.csv'), R);
write_sample_index(fullfile(AF_DIR, 'dataset_index.csv'), samp_info);
fig_summary(R, fullfile(FIG_DIR, 'autofocus_summary.png'));
fig_case(R(cellfun(@(x) strcmp(x, 'E2'), {R.id})), CASE_DIR, fullfile(FIG_DIR, 'autofocus_case_quad.png'));
fig_case(R(cellfun(@(x) strcmp(x, 'E5'), {R.id})), CASE_DIR, fullfile(FIG_DIR, 'autofocus_case_multitone.png'));
fig_case(R(cellfun(@(x) strcmp(x, 'E3'), {R.id})), CASE_DIR, fullfile(FIG_DIR, 'autofocus_case_sin.png'));
fig_entropy_curves(R, fullfile(FIG_DIR, 'autofocus_entropy_curves.png'));
fprintf('  CSV: %s\n', fullfile(MET_DIR, 'autofocus_metrics.csv'));
fprintf('  索引: %s\n', fullfile(AF_DIR, 'dataset_index.csv'));

%% 6. 验收（规划方案 七.1 自聚焦条目）
fprintf('\n=== 验收结论（规划方案 七.1 · 自聚焦）===\n');
main_ids  = {'E1','E2','E3','E4','E5','E6','E7','E8'};
zero_ids  = {'E0'};
unseen_ids= {'E9'};

get = @(id) R(strcmp({R.id}, id));
n_main = 0; n_meaf = 0; n_pga = 0;
for k = 1:numel(main_ids)
    rr = get(main_ids{k});
    if isempty(rr); continue; end
    n_main = n_main + 1;
    if rr.nrmse_meaf < rr.nrmse_def; n_meaf = n_meaf + 1; end
    if rr.nrmse_pga  < rr.nrmse_def; n_pga  = n_pga  + 1; end
end

fprintf('非零误差案例数 = %d（要求 >= 6）\n', n_main);
fprintf('NRMSE 优于散焦输入的案例数：MEAF %d/%d，PGA %d/%d（要求至少 5）\n', ...
    n_meaf, n_main, n_pga, n_main);

ok_main = (n_main >= 6) && (n_meaf >= 5 || n_pga >= 5);
fprintf('主实验（>=6 案例且 >=5 改善）: %s\n', cond_str(ok_main));

ok_zero = true;
for k = 1:numel(zero_ids)
    rr = get(zero_ids{k});
    if isempty(rr); continue; end
    % 零误差：算法不应显著劣化（NRMSE 增量 < 0.02，熵增幅 < 0.02 bit）
    d_nrmse = rr.nrmse_meaf - rr.nrmse_def;
    d_ent   = rr.entropy_meaf - rr.entropy_def;
    good = (d_nrmse < 0.02) && (d_ent < 0.02);
    fprintf('[零误差 %s] NRMSE %.4f -> %.4f (Δ%+.4f)，熵 %.4f -> %.4f (Δ%+.4f) : %s\n', ...
        rr.id, rr.nrmse_def, rr.nrmse_meaf, d_nrmse, ...
        rr.entropy_def, rr.entropy_meaf, d_ent, cond_str(good));
    ok_zero = ok_zero && good;
end
fprintf('零误差鲁棒性: %s\n', cond_str(ok_zero));

ok_unseen = true;
for k = 1:numel(unseen_ids)
    rr = get(unseen_ids{k});
    if isempty(rr); continue; end
    good = rr.nrmse_meaf < rr.nrmse_def || rr.nrmse_pga < rr.nrmse_def;
    fprintf('[未见误差 %s] %s | NRMSE 散焦=%.4f MEAF=%.4f PGA=%.4f : %s\n', ...
        rr.id, rr.desc, rr.nrmse_def, rr.nrmse_meaf, rr.nrmse_pga, cond_str(good));
    ok_unseen = ok_unseen && good;
end
fprintf('未见误差泛化: %s\n', cond_str(ok_unseen));

t_meaf_all = mean([R.time_meaf_s]);
t_pga_all  = mean([R.time_pga_s]);
fprintf('平均耗时（%dx%d 单块）：MEAF %.2f s / PGA %.2f s\n', Na, Na, t_meaf_all, t_pga_all);

rep = struct();
rep.generated_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
rep.block_size = [BLK BLK];
rep.n_nonzero_cases = n_main;
rep.n_nrmse_improved_meaf = n_meaf;
rep.n_nrmse_improved_pga = n_pga;
rep.acceptance_main = ok_main;
rep.acceptance_zero = ok_zero;
rep.acceptance_unseen = ok_unseen;
rep.mean_time_meaf_s = t_meaf_all;
rep.mean_time_pga_s = t_pga_all;
rep.cases = rmfield(R, 'entropy_hist_pga');
fid = fopen(fullfile(AF_DIR, 'autofocus_week3_report.json'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s', jsonencode(rep, 'PrettyPrint', true));
fclose(fid);

fprintf('\n全部验收: %s\n', cond_str(ok_main && ok_zero && ok_unseen));
fprintf('报告: %s\n', fullfile(AF_DIR, 'autofocus_week3_report.json'));


%% ============================ Local Functions ============================

function bases = build_bases(sim_mat, slc_mat, focus_mat, blk)
% 构造三张 512x512 参考聚焦像：仿真多点目标 / 真实 SM SLC / 真实 Level-0 聚焦。
    bases = struct();

    % --- B1 仿真多点目标：裁到含 (513,701) 目标的 512x512 窗口，目标落在中心附近
    S = load(sim_mat, 'img_m');
    img_m = S.img_m;
    r0 = 513 - floor(blk/2); c0 = 701 - floor(blk/2);
    r0 = max(1, min(size(img_m,1)-blk+1, r0));
    c0 = max(1, min(size(img_m,2)-blk+1, c0));
    B1 = double(img_m(r0:r0+blk-1, c0:c0+blk-1));
    bases.B1 = struct('img', B1, 'src', sprintf('sim_multi[%d:%d,%d:%d]', r0, r0+blk-1, c0, c0+blk-1), ...
        'spikiness', spikiness(B1), 'kind', 'simulated');

    % --- B2 真实 SM SLC：候选 512 块（步长 256）中取尖峰度最高者（自聚焦需要强散射体）
    S = load(slc_mat, 'slc');
    slc = double(S.slc);
    best = []; best_sp = -Inf;
    step = 256;
    for rr = 1:step:(size(slc,1)-blk+1)
        for cc = 1:step:(size(slc,2)-blk+1)
            w = slc(rr:rr+blk-1, cc:cc+blk-1);
            sp = spikiness(w);
            if sp > best_sp; best_sp = sp; best = [rr cc]; end
        end
    end
    if isempty(best); best = [1 1]; end
    B2 = slc(best(1):best(1)+blk-1, best(2):best(2)+blk-1);
    bases.B2 = struct('img', B2, 'src', sprintf('slc_VV[%d:%d,%d:%d]', best(1), best(1)+blk-1, best(2), best(2)+blk-1), ...
        'spikiness', spikiness(B2), 'kind', 'real_slc');

    % --- B3 真实 Level-0 聚焦结果（v7.3 大文件，分段局部读）
    mf = matfile(focus_mat);
    [Nr, Nc] = size(mf, 'img_real');
    r_cands = 1:4096:max(1, Nr-blk+1);
    c_cands = unique([1, max(1, floor((Nc-blk)/2)), max(1, Nc-blk+1)]);
    best = []; best_sp = -Inf;
    for rr = r_cands
        for cc = c_cands
            w = double(mf.img_real(rr:rr+blk-1, cc:cc+blk-1));
            sp = spikiness(w);
            if sp > best_sp; best_sp = sp; best = [rr cc]; end
        end
    end
    if isempty(best); best = [1 1]; end
    B3 = double(mf.img_real(best(1):best(1)+blk-1, best(2):best(2)+blk-1));
    bases.B3 = struct('img', B3, 'src', sprintf('level0_focused[%d:%d,%d:%d]', best(1), best(1)+blk-1, best(2), best(2)+blk-1), ...
        'spikiness', spikiness(B3), 'kind', 'real_level0');
end

function s = spikiness(img)
% 尖峰度 = 幅度 99.9 分位 / 中位数，用于挑选含强散射体的块。
    a = abs(img(:));
    s = my_prctile(a, 99.9) / (my_prctile(a, 50) + 1e-12);
end

function cases = define_cases()
% 8 个独立非零误差 + 1 个零误差 + 1 个未见误差。
% 误差类型：quad 二次 / sin 单频正弦 / multitone 多频正弦 / rand 带限随机 / mixed 混合
    cases = struct('id', {}, 'base', {}, 'type', {}, 'desc', {}, 'seed', {}, 'p', {});
    add = @(id, base, type, desc, seed, p) struct('id', id, 'base', base, 'type', type, ...
        'desc', desc, 'seed', seed, 'p', p);

    cases(end+1) = add('E1', 'B1', 'quad',     '二次相位 a2=6 rad（Ka 失配·弱）',      11, struct('a2', 6));
    cases(end+1) = add('E2', 'B1', 'quad',     '二次相位 a2=14 rad（Ka 失配·强）',     12, struct('a2', 14));
    cases(end+1) = add('E3', 'B1', 'sin',      '单频正弦 A=1.2 rad·1.5 周期',          13, struct('A', 1.2, 'n', 1.5));
    cases(end+1) = add('E4', 'B2', 'quad',     '二次相位 a2=8 rad（真实 SLC）',        14, struct('a2', 8));
    cases(end+1) = add('E5', 'B2', 'multitone', '双频正弦 A=1.0 rad·1+3.5 周期',       15, struct('A', 1.0, 'n1', 1.0, 'n2', 3.5));
    cases(end+1) = add('E6', 'B2', 'rand',     '带限随机 rms=0.7 rad·平滑宽 8 bin',    16, struct('rms', 0.7, 'w', 8));
    cases(end+1) = add('E7', 'B3', 'sin',      '单频正弦 A=1.5 rad·2.5 周期（L0）',    17, struct('A', 1.5, 'n', 2.5));
    cases(end+1) = add('E8', 'B3', 'rand',     '高频随机 rms=1.0 rad·平滑宽 2 bin',    18, struct('rms', 1.0, 'w', 2));
    cases(end+1) = add('E0', 'B1', 'zero',     '零误差（无注入）',                     19, struct());
    cases(end+1) = add('E9', 'B2', 'mixed',    '未见：二次5+正弦0.8+随机0.4 混合',      20, ...
        struct('a2', 5, 'A', 0.8, 'n', 1.8, 'rms', 0.4, 'w', 5));
end

function u = azimuth_u(Na)
% 方位频率归一化坐标（与 fft 的原始排序一致，DC 在 index 1），取值 [-1, 1]。
    u = [(0:Na/2-1), (-Na/2:-1)] / (Na/2);
end

function phi = make_error(u, cs, seed)
% 按案例定义生成 1 x Na 的已知相位误差（弧度）。
    Na = numel(u);
    switch cs.type
        case 'zero'
            phi = zeros(1, Na);
        case 'quad'
            phi = cs.p.a2 * u.^2;
        case 'sin'
            phi = cs.p.A * sin(2 * pi * cs.p.n * (u + 1) / 2);
        case 'multitone'
            phi = cs.p.A * (0.6 * sin(2 * pi * cs.p.n1 * (u + 1) / 2) + ...
                            0.4 * sin(2 * pi * cs.p.n2 * (u + 1) / 2));
        case 'rand'
            phi = bandlimited_random(u, cs.p.w, cs.p.rms, seed);
        case 'mixed'
            phi = cs.p.a2 * u.^2 ...
                + cs.p.A * sin(2 * pi * cs.p.n * (u + 1) / 2) ...
                + bandlimited_random(u, cs.p.w, cs.p.rms, seed);
        otherwise
            error('未知误差类型: %s', cs.type);
    end
end

function phi = bandlimited_random(u, width_bins, rms_target, seed)
% 带限随机相位：白噪声经循环高斯平滑（宽度 width_bins）后归一到指定 RMS。
    Na = numel(u);
    rs = rng;                       % 保存并恢复全局随机流，保证可复现且不串扰
    rng(seed, 'twister');
    g = randn(1, Na);
    n = max(1, round(width_bins));
    t = -(n-1):(n-1);
    ker = exp(-(t.^2) / (2 * (max(n,2)/2)^2));
    ker = ker / sum(ker);
    gc = conv(g, ker, 'same');      % 线性卷积取 same，边缘影响可忽略
    rng(rs);
    gc = gc - mean(gc);
    phi = rms_target * gc / (sqrt(mean(gc.^2)) + 1e-12);
end

function out = apply_phase(img, phi)
% 在方位频谱（第 2 维）上乘 exp(1j*phi)。传 +phi 为「注入」，传 -phi 为「补偿」。
    phi = phi(:).';
    S = fft(img, [], 2);
    S = S .* exp(1j * phi);
    out = ifft(S, [], 2);
end

function B = meaf_basis(u)
% 最小熵自聚焦的相位基：
%   偶次多项式 u^2,u^4,u^6 —— 描述聚焦参数失配（二次/高阶散焦）
%   余弦/正弦各 4 阶 —— 描述单频、多频振荡型误差（必须同时含 cos 与 sin，否则
%   奇对称误差无法表达；实测只给 cos 时正弦误差几乎修不回来）
    u = u(:);
    cols = {u.^2, u.^4, u.^6};
    for m = 1:4
        cols{end+1} = cos(2*pi*m*(u+1)/2);   %#ok<AGROW>
        cols{end+1} = sin(2*pi*m*(u+1)/2);   %#ok<AGROW>
    end
    B = [cols{:}];
end

function [img_c, phi_hat, info] = meaf_autofocus(img, u, cfg)
% 最小熵自聚焦：以图像熵为目标，在相位基上做 Nelder-Mead 无约束极小化。
    B = meaf_basis(u);                  % Na x K
    K = size(B, 2);
    scale = [3 2 2 ones(1, K-3) * 1.5]; % 各基函数系数的量级（rad）
    obj = @(x) image_entropy(apply_phase(img, -bfun(B, scale, x)));

    % 1) 粗扫二次项（散焦主因），给 fminsearch 一个好初值
    a2_grid = linspace(-cfg.coarse_max, cfg.coarse_max, cfg.coarse_n);
    e_grid = zeros(size(a2_grid));
    for i = 1:numel(a2_grid)
        c = zeros(K, 1); c(1) = a2_grid(i);
        e_grid(i) = image_entropy(apply_phase(img, -(B * c)));
    end
    [~, im] = min(e_grid);
    x0 = 1e-3 * ones(K, 1);
    x0(1) = a2_grid(im) / scale(1);

    % 2) Nelder-Mead 精调
    opts = optimset('MaxFunEvals', cfg.maxfev, 'MaxIter', cfg.maxfev, ...
                    'Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-4);
    x = fminsearch(obj, x0, opts);
    c = scale(:) .* x(:);
    phi_hat = (B * c).';
    img_c = apply_phase(img, -phi_hat);
    info = struct('coeff', c, 'coarse_best_a2', a2_grid(im));
end

function p = bfun(B, scale, x)
% 相位基线性组合（便于匿名函数调用）。
    p = B * (scale(:) .* x(:));
end

function W = wrap_hann(Na, w)
% 以 index 1 为峰值、环绕的升余弦窗（宽 w 采样）。
    w = max(3, round(w));
    m = -floor(w/2):floor(w/2);
    wts = 0.5 + 0.5 * cos(2 * pi * m / w);
    W = zeros(1, Na);
    W(mod(m, Na) + 1) = wts;
end

function [img_c, phi_hat, ehist, n_acc] = pga_autofocus(img, cfg)
% 相位梯度自聚焦（PGA）。谱域相位只在有效频带内 unwrap/去线性；加熵不劣化守卫。
    [Nr, Na] = size(img);
    e = sum(abs(img).^2, 2);
    thr = my_prctile(e, 100 * (1 - cfg.bin_frac));
    bins = find(e >= thr);
    if numel(bins) < cfg.min_bins
        [~, ord] = sort(e, 'descend');
        bins = ord(1:min(cfg.min_bins, Nr));
    end
    bins = bins(:).';

    img_c = img;
    phi_hat = zeros(1, Na);
    h_now = image_entropy(img_c);
    nIt = numel(cfg.win_frac);
    ehist = zeros(1, nIt);
    n_acc = 0;

    for it = 1:nIt
        W = wrap_hann(Na, Na * cfg.win_frac(it));
        X = zeros(numel(bins), Na);
        for j = 1:numel(bins)
            row = img_c(bins(j), :);
            [~, pk] = max(abs(row));
            row = circshift(row, [0, 1 - pk]);       % 亮点圆移到 index 1（不可改成中心）
            X(j, :) = fft(row .* W);
        end
        Xs = fftshift(X, 2);
        magc = sum(abs(Xs), 1);
        valid = magc > cfg.band_thr * max(magc);
        kk = find(valid);
        if numel(kk) < 8; break; end

        a = unwrap(angle(sum(Xs(:, valid), 1)));      % 只在带内 unwrap
        a = a - detrend_line(kk(:), a(:));            % 只在带内去线性
        ph = zeros(1, Na);
        ph(kk) = a;
        ph(1:kk(1)-1) = a(1);                         % 带外保持常数，避免边界跳变
        ph(kk(end)+1:end) = a(end);

        cand = apply_phase(img_c, -ifftshift(ph));
        h_cand = image_entropy(cand);
        if h_cand < h_now                              % 熵不劣化守卫
            img_c = cand;
            phi_hat = phi_hat + ifftshift(ph);
            h_now = h_cand;
            n_acc = n_acc + 1;
            ehist(it) = h_cand;
        else
            ehist(it) = h_now;
            break;
        end
    end
    phi_hat = phi_hat - remove_linear(phi_hat);
end

function y = detrend_line(x, v)
% 在给定样本点 x 上用一次多项式拟合 v，返回拟合值（行向量）。
    A = [ones(numel(x), 1) x(:)];
    c = A \ v(:);
    y = (A * c).';
end

function y = remove_linear(x)
% 去掉常数项与线性项（最小二乘拟合后扣除）。
    x = x(:);
    n = numel(x);
    k = (1:n).';
    A = [ones(n, 1) k];
    c = A \ x;
    y = x - A * c;
    y = y(:).';
end

function r = phase_resid_rms(phi_hat, phi_true)
% 估计相位与真值的残差 RMS（先扣除常数与线性项，二者不影响聚焦）。
    d = phi_hat(:) - phi_true(:);
    d = remove_linear(d);
    r = sqrt(mean(d.^2));
end

function e = image_entropy(image)
% 归一化香农熵（基于幅度平方分布），单位 bit。
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

function v = nrmse_amp(ref, test)
% 幅度 NRMSE = || |test| - |ref| ||_F / || |ref| ||_F。
    a = abs(ref); b = abs(test);
    v = norm(b(:) - a(:)) / (norm(a(:)) + 1e-12);
end

function [w_ref, w_def, w_meaf, w_pga] = azimuth_width3db(ref, def, img_meaf, img_pga)
% 过参考图峰值点的方位向 -3dB 主瓣宽度（采样数），三点线性插值。
    [~, idx] = max(abs(ref(:)));
    [pr, ~] = ind2sub(size(ref), idx);
    w_ref  = width3db(abs(ref(pr, :)));
    w_def  = width3db(abs(def(pr, :)));
    w_meaf = width3db(abs(img_meaf(pr, :)));
    w_pga  = width3db(abs(img_pga(pr, :)));
end

function w = width3db(prof)
% 单条剖面 -3dB 宽度（采样数）。
    prof = double(prof(:)).';
    n = numel(prof);
    [pk_val, pk] = max(prof);
    if pk_val <= 0; w = NaN; return; end
    half = pk_val / sqrt(2);
    i = pk;
    while i > 1 && prof(i) > half; i = i - 1; end
    if i >= 1 && i < n && prof(i+1) ~= prof(i)
        left3 = i + (half - prof(i)) / (prof(i+1) - prof(i));
    else
        left3 = i;
    end
    i = pk;
    while i < n && prof(i) > half; i = i + 1; end
    if i > 1 && i <= n && prof(i) ~= prof(i-1)
        right3 = i - 1 + (half - prof(i-1)) / (prof(i) - prof(i-1));
    else
        right3 = i;
    end
    w = right3 - left3;
end

function info = make_training_samples(slc_mat, out_dir, blk, n_win, n_err)
% 从真实 SM SLC 抽 blk x blk 块，每块注入 n_err 种随机/正弦误差，保存网络训练样本对。
    S = load(slc_mat, 'slc');
    slc = double(S.slc);
    info = struct([]);
    rng(20261107, 'twister');
    used = 0;
    tries = 0;
    while used < n_win && tries < 200
        tries = tries + 1;
        r0 = randi([1, size(slc,1) - blk + 1]);
        c0 = randi([1, size(slc,2) - blk + 1]);
        ref = slc(r0:r0+blk-1, c0:c0+blk-1);
        ref = ref - mean(ref(:));                 % 去直流，避免零频主导熵
        % 只保留含足够强散射体的窗口
        if spikiness(ref) < 8; continue; end
        used = used + 1;
        u = azimuth_u(blk);
        for e = 1:n_err
            switch mod(e, 2)
                case 1
                    cs = struct('type', 'quad', 'p', struct('a2', 4 + 10 * rand()));
                case 0
                    cs = struct('type', 'rand', 'p', struct('rms', 0.4 + 0.8 * rand(), 'w', 2 + 10 * rand()));
            end
            phi = make_error(u, cs, 1000 + used * 10 + e);
            def = apply_phase(ref, -phi);
            fn = sprintf('af_pair_w%02d_e%d.mat', used, e);
            defocus = single(def); reference = single(ref);
            save(fullfile(out_dir, fn), 'defocus', 'reference', 'phi', 'cs', 'r0', 'c0', '-v7.3');
            s = struct('file', fn, 'type', cs.type, 'phase_rms', sqrt(mean(phi.^2)), ...
                       'row0', r0, 'col0', c0, 'size', blk);
            if isempty(info); info = s; else; info(end+1) = s; end %#ok<SAGROW>
        end
    end
end

function write_metrics_csv(path, R)
% 逐案例指标表。
    fid = fopen(path, 'w', 'n', 'UTF-8');
    fprintf(fid, 'case,base,error_type,phase_rms_rad,entropy_ref,entropy_defocus,entropy_meaf,entropy_pga,');
    fprintf(fid, 'contrast_ref,contrast_defocus,contrast_meaf,contrast_pga,');
    fprintf(fid, 'nrmse_defocus,nrmse_meaf,nrmse_pga,');
    fprintf(fid, 'az_width3db_ref,az_width3db_defocus,az_width3db_meaf,az_width3db_pga,');
    fprintf(fid, 'time_meaf_s,time_pga_s,phi_resid_meaf_rad,phi_resid_pga_rad,desc\n');
    for k = 1:numel(R)
        r = R(k);
        fprintf(fid, '%s,%s,%s,%.4f,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%.4f,%.4f,', ...
            r.id, r.base, r.type, r.phi_rms_inj, r.entropy_ref, r.entropy_def, ...
            r.entropy_meaf, r.entropy_pga, r.contrast_ref, r.contrast_def, ...
            r.contrast_meaf, r.contrast_pga);
        fprintf(fid, '%.6f,%.6f,%.6f,%.4f,%.4f,%.4f,%.4f,%.3f,%.3f,%.4f,%.4f,%s\n', ...
            r.nrmse_def, r.nrmse_meaf, r.nrmse_pga, r.width_ref, r.width_def, ...
            r.width_meaf, r.width_pga, r.time_meaf_s, r.time_pga_s, ...
            r.phi_resid_meaf, r.phi_resid_pga, r.desc);
    end
    fclose(fid);
end

function write_sample_index(path, info)
% 样本索引 CSV（供第 4 周训练脚本按场景/误差类型读取）。
    fid = fopen(path, 'w', 'n', 'UTF-8');
    fprintf(fid, 'file,error_type,phase_rms_rad,row0,col0,size\n');
    for k = 1:numel(info)
        s = info(k);
        fprintf(fid, '%s,%s,%.4f,%d,%d,%d\n', s.file, s.type, s.phase_rms, s.row0, s.col0, s.size);
    end
    fclose(fid);
end

function fig_summary(R, path)
% 汇总图：左=各案例熵（参考/散焦/MEAF/PGA），右=各案例幅度 NRMSE。
    n = numel(R);
    ids = {R.id};
    E = [[R.entropy_ref]; [R.entropy_def]; [R.entropy_meaf]; [R.entropy_pga]].';
    N = [[R.nrmse_def]; [R.nrmse_meaf]; [R.nrmse_pga]].';
    h = figure('Visible', 'off', 'Position', [100 100 1280 480]);
    subplot(1, 2, 1);
    bar(E); grid on;
    set(gca, 'XTick', 1:n, 'XTickLabel', ids);
    legend({'reference', 'defocused', 'MEAF', 'PGA'}, 'Location', 'northwest');
    ylabel('image entropy (bit)'); title('Entropy per case'); xlabel('case');
    subplot(1, 2, 2);
    bar(N); grid on;
    set(gca, 'XTick', 1:n, 'XTickLabel', ids);
    legend({'defocused', 'MEAF', 'PGA'}, 'Location', 'northwest');
    ylabel('amplitude NRMSE'); title('NRMSE vs reference'); xlabel('case');
    print(path, '-dpng', '-r120'); close(h);
end

function fig_case(R, case_dir, path)
% 单案例 2x2 对比图：参考 / 散焦 / MEAF / PGA（对数幅度，统一显示范围）。
    if isempty(R); return; end
    S = load(fullfile(case_dir, sprintf('af_case_%s.mat', R.id)));
    imgs = {S.ref, S.def, S.img_meaf, S.img_pga};
    tits = {'reference (focused)', 'defocused input', 'MEAF corrected', 'PGA corrected'};
    db = cell(1, 4);
    allv = [];
    for k = 1:4
        a = abs(imgs{k}); a(a <= 0) = 1e-12;
        db{k} = 20 * log10(a / max(a(:)));
        allv = [allv; db{k}(:)];
    end
    vmin = my_prctile(allv, 5); vmax = 0;
    h = figure('Visible', 'off', 'Position', [100 100 1000 820]);
    for k = 1:4
        subplot(2, 2, k);
        imagesc(db{k}, [vmin vmax]); colormap('gray'); axis image;
        title(sprintf('%s\nNRMSE=%.4f  H=%.4f', tits{k}, nrmse_pick(R, k), entropy_pick(R, k)));
        if k >= 3; xlabel('azimuth (samples)'); end
        if mod(k, 2) == 1; ylabel('range (samples)'); end
    end
    if ~isempty(which('sgtitle'))
        sgtitle(sprintf('%s  %s', R.id, R.desc), 'Interpreter', 'none');
    end
    print(path, '-dpng', '-r110'); close(h);
end

function v = nrmse_pick(R, k)
    switch k
        case 1; v = 0;
        case 2; v = R.nrmse_def;
        case 3; v = R.nrmse_meaf;
        case 4; v = R.nrmse_pga;
    end
end

function v = entropy_pick(R, k)
    switch k
        case 1; v = R.entropy_ref;
        case 2; v = R.entropy_def;
        case 3; v = R.entropy_meaf;
        case 4; v = R.entropy_pga;
    end
end

function fig_entropy_curves(R, path)
% PGA 逐迭代熵曲线（体现收敛性）。
    h = figure('Visible', 'off', 'Position', [100 100 900 460]);
    hold on; grid on;
    for k = 1:numel(R)
        eh = R(k).entropy_hist_pga;
        plot(0:numel(eh)-1, eh, '-o', 'MarkerSize', 3);
    end
    xlabel('PGA iteration'); ylabel('image entropy (bit)');
    legend({R.id}, 'Location', 'eastoutside');
    title('PGA entropy convergence');
    print(path, '-dpng', '-r120'); close(h);
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
