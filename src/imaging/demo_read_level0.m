% demo_read_level0.m — read_level0 用法示例（A组 第2周 Level-0 读取）
%
% 运行方式：
%   1) 在 MATLAB 中直接运行本脚本；
%   2) 命令行： matlab -batch "run('D:/study/SAR/src/imaging/demo_read_level0.m')"
%
% 说明：RawRoot 指向助教提供的 Level-0 RAW 产品根目录（含 *SAFE）。
%       本产品为 SDV（VV/VH），务必指定 Polarization='VV'（或 'VH'）。

rawRoot = 'D:/study/SAR/data/raw/S1A_S2_RAW__0SDV_20250303T140022_20250303T140047_058141_072E7E_7DD0';

% 1) 小块验证：解码 1500 行 -> IQ 校正 + SWST 对齐 -> 校正块/参数 JSON/3 个 ROI
%    Plot=false 便于无显示环境（批处理/服务器）运行；有桌面时可改 true 生成 ROI 预览
read_level0('RawRoot', rawRoot, 'Polarization', 'VV', 'NumPulses', 1500, 'Plot', false);

% 2) 整景磁盘模式（约 10GB MAT，耗时；需要时取消下一行注释）
% read_level0('RawRoot', rawRoot, 'Polarization', 'VV', 'Full', true);
