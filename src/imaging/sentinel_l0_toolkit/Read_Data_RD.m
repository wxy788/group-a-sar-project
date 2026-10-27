clear;clc;close all;
%% -------- 数据读取 --------
%% 参数
Na                              = 16384;
startPulse                      = 60644/2-Na/2;
Nfast                           = 16384;
fc                              = 5.405e9;
C                               = 299792458;

%% 方位选择
inputFile = ['D:\study\SAR\data\S1A_S2_RAW__0SDV_20250303T140022_20250303T140047_058141_072E7E_7DD0\S1A_S2_RAW__0SDV_20250303T140022_20250303T140047_058141_072E7E_7DD0.SAFE\s1a-s2-raw-s-vh-20250303t140022-20250303t140047-058141-072e7e.dat'];
[echo, info] = read_s1_raw(inputFile, Na, ...
    'StartPulse', startPulse, ...
    'CarrierFrequencyHz', fc, ...
    'Polarization', 'VH', ...
    'Plot', false);

[echo, info]                    = IQ_correction(echo, info);   % 先对原始有效样本进行 IQ 校正
[echo, info]                    = echo_SWST(echo, info);       % 再对齐各脉冲的接收窗口
[echo, info]                    = Shift_Delay_correction(echo, info, inputFile); % 漂移与内部时延补偿
info                            = doppler_geo_parameter(info); % XML 多普勒、轨道与几何参数

Nr                              = size(echo,2);
r_start                         = Nr/2-Nfast/2+1;
r_end                           = Nr/2+Nfast/2;

fc                              = info.radar.fc;
lambda                          = info.radar.lambda;
Fr                              = info.radar.range_sample_rate_hz;
PRF                             = info.radar.prf_hz;

tmp                             = info.packet_metadata.values;
idx                             = tmp.Signal_Type == 0 & tmp.selected_row_zero_based >= 0;
Tp                              = unique(tmp.Tx_Pulse_Length(idx));
Kr                              = unique(tmp.Tx_Ramp_Rate(idx));
Br                              = unique(abs(tmp.Tx_Ramp_Rate(idx)).*tmp.Tx_Pulse_Length(idx));
tr                              = info.radar.swst_reference_s+info.radar.decimation_transient_delay_s+(0:Nr-1)/Fr;
rSlant                          = info.radar.slant_range_m;
ta                              = info.radar.pulse_time_relative_s;
ta                              = ta-ta(end/2-1);
S_ECEF                          = info.orbit.pulse_position_ecef_m;
vS_ECEF                         = info.orbit.pulse_velocity_ecef_mps;
qS                              = info.attitude.quaternion;
fa                              = (-Na/2:Na/2-1)*PRF/Na;
Np                              = ceil(Tp*Fr);

fc_start                        = tmp.Tx_Pulse_Start_Frequency(idx);

figure;
imagesc(tr,ta,abs(echo));
xlabel('快时间');ylabel('慢时间');title('回波矩阵');

%% 距离选择
DATA                            = echo.';
DATA                            = DATA(r_start:r_end,:);
tr                              = tr(r_start:r_end).';
rSlant                          = rSlant(r_start:r_end).';
Nr                              = size(DATA,1);
clear echo;

%% ------------------ Main_RADARSAT_1_Strimap_RDA.m -----------------------
% 参数
Vr                              = info.geometry.center_equivalent_velocity_mps;
dR                              = info.radar.range_sample_spacing_m;
Tr0                             = info.radar.slant_range_time_s(1);
R0                              = info.radar.slant_range_m(1);
tnrn                            = (-Nr/2:Nr/2-1)/Fr;
fr                              = (-Nr/2:Nr/2-1)/Nr*Fr;
Ka                              = info.doppler.center_azimuth_fm_rate_hz_s;
Fdc                             = info.doppler.center_centroid_hz;
SquintAngle                     = asin(lambda*Fdc/2/Vr);
AzimuthBeamWidth                = lambda*PRF/2/Vr/cos(SquintAngle);
Rs                              = info.radar.center_slant_range_m;
dA                              = Vr/PRF;
xCenter                         = Rs*Fdc*lambda/2/Vr;
yCenter                         = sqrt(Rs^2-xCenter^2);
La                              = Rs*AzimuthBeamWidth/cos(SquintAngle);
TLa                             = La/Vr;
Ba                              = abs(TLa*Ka);
NLa                             = floor(La/dA);
ta                              = ta.';

% 距离多普勒成像算法处理 精确RCMC+方式2(精确方式)SRC
% -------------------------------------------------------------------------
% DATA      : 转置后的原始数据(列为距离,行为方位)
% SARPara   : SAR系统参数
% ImagePara : 成像参数
%
% DATA      : 成像数据
% -------------------------------------------------------------------------
% 2022-12-02  Junxu Wang  创建
% 2022-12-02  Junxu Wang  完成调试
%                         为了补零不混叠,作FFT的次数比理论的上多
% 2022-12-04  Junxu Wang  方位匹配滤波器与距离相关
% 2022-12-07  Junxu Wang  SRC滤波器改为固定距离
% 2022-12-24  Junxu Wang  修正二次距离压缩方式
% -------------------------------------------------------------------------
% TODO:考虑校正斜视角引起的方位偏移. (待确认)

Rref                            = Rs*cos(SquintAngle);
NfftR                           = 2^nextpow2(Np + Nr - 1);
fa                              = fa + Fdc;

% 脉冲压缩
NHf                             = (-NfftR / 2 : NfftR / 2 - 1)' / NfftR * Fr;
Nzero                           = fix((Fr - Br) / Fr * NfftR / 2) * 2;
% WinR                            = padarray(kaiser(NfftR - Nzero, 2.5), [Nzero/2 0], 0);
WinR                            = padarray(ones(NfftR - Nzero, 1), [Nzero/2 0], 0);
Hf                              = WinR .* exp(1j * pi * NHf.^2 ./ Kr);

% % 为了处理大数据,采用循环而不用矩阵运算减低对内存的要求
% for aa = 1:Na
%     if mod(aa, 128) == 0
%         fprintf('RC Processing: %d / %d\n', aa, Na);
%     end
% 
%     temp   = DATA(:,aa) * exp(-1j*2*pi*Fdc*ta(aa));     % 方位多普勒校正
%     temp   = padarray(temp, [(NfftR - Nr)/2 0], 0);     % 回波数据两边补零
%     temp   = fftshift(ifft(...                          % 脉冲压缩
%                          fft(ifftshift(temp, 1)) .* ifftshift(Hf, 1)), 1);
%     DATA(:,aa) = temp(NfftR / 2 - Nr / 2 + 1: ...          % 舍去补零区
%                       NfftR / 2 + Nr / 2);
% end

temp                            = DATA .* exp(-1j*2*pi*Fdc*ta);
temp                            = padarray(temp, [(NfftR - Nr)/2 0], 0);
temp                            = fftshift(ifft(fft(ifftshift(temp, 1)) .* ifftshift(Hf, 1)), 1);
DATA                            = temp(NfftR / 2 - Nr / 2 + 1 : NfftR / 2 + Nr / 2, :);
clear temp Hf;

% figure; imagesc(abs(DATA));
% title('距离压缩结果');

% figure; imagesc(abs(fftshift(fft(DATA,[],2),2)));
% title('距离多普勒结果');

% 二次距离压缩
% -------------------------------------------------------------------------
DATA                            = fftshift(fft2(ifftshift(DATA)));
D                               = sqrt(1-(lambda*fa/2/Vr).^2); % 徙动因子[式(6.24)]
Ksrc                            = 2*Vr^2*fc^3*D.^3./(C*Rref*fa.^2); % SRC滤波器调频率[式(6.22)]
fr                              = fr.';

Hsrc                            = exp(-1j*pi*fr.^2./Ksrc);
DATA                            = DATA.*Hsrc;

DATA                            = fftshift(ifft(ifftshift(DATA, 1), [], 1), 1);
clear Hsrc;

% figure; imagesc(abs(DATA));
% title('二次距离压缩结果');

% 距离徙动校正
% 徙动校正(精确实现)
% -------------------------------------------------------------------------
% Rrefs  = (Rs + Tnrn*C/2) * cos(SquintAngle);            % 每个距离单元的Rref 
% RCM    = Rrefs.'./D - (Rs + Tnrn*C/2)';                 % 距离徙动量[式(6.25)],对齐场景中心

% 2026/10/06 york 校正到最近斜距R_0
Rgrid                           = Rs + tnrn*C/2;
Rrefs                           = Rgrid;
RCM                             = Rrefs.'./D - Rgrid.';

N_RCM                           = ceil(RCM/dR);   % 距离徙动单元(整数部分)
F_RCM                           = N_RCM - RCM/dR; % 距离徙动单元(小数部分)
NKel                            = 8;              % 插值核点数

DATA_RCMC                       = zeros(Nr, Na);

for aa       = 1:Na
    if mod(aa,128) == 0
        fprintf('RCMC Processing: %d / %d\n', aa, Na);
    end

    for rr   = 1:Nr
        % sinc插值核
        SincKel = sinc(F_RCM(rr,aa)-NKel/2:1:F_RCM(rr,aa)+NKel/2-1);
        SincKel = SincKel/sum(SincKel);                 % 插值核归一化
        kk      = rr + N_RCM(rr,aa);                    % 原矩阵中距离徙动的位置
        kk      = rem(kk, Nr);                          % 整数部分超过Nr的话取余

        % 插值核的取值位置
        % -----------------------------------------------------------------
        if     (kk - NKel/2)     > Nr                   % 全右溢出
            nn  = (kk - NKel/2 - Nr):1:(kk + NKel/2 - 1 - Nr);
        elseif (kk + NKel/2 - 1) > Nr                   % 部分右溢
            nn1 = (kk - NKel/2):1:Nr;
            nn2 = 1:1:(kk + NKel/2 -1 - Nr);
            nn  = [nn1, nn2];
        elseif (kk + NKel/2 - 1) < 1                    % 全左溢出
            nn  = (kk - NKel/2 + Nr):1:(kk + NKel/2 - 1 +Nr);
        elseif (kk - NKel/2)     < 1                    % 部分左溢
            nn1 = (kk - NKel/2 + Nr):1:Nr;
            nn2 = 1:1:(kk + NKel/2 - 1);
            nn  = [nn1,nn2];
        else                                            % 无溢出
            nn  = (kk - NKel/2):1:(kk + NKel/2 -1);
        end

        DATA_RCMC(rr,aa) = SincKel*DATA(nn,aa);         % 计算插值
    end
end

DATA  = DATA_RCMC; clear DATA_RCMC

clear N_RCM F_RCM RCM
% figure; imagesc(abs(DATA));
% title('距离徙动校正后RD域');

% 方位压缩
% -------------------------------------------------------------------------
NfftA                           = 2^nextpow2(NLa + Na - 1); % 方位向FFT点数
Fnan                            = Fdc*1 + (-NfftA/2:NfftA/2-1)/NfftA * PRF;
fbase                           = (-NfftA/2:NfftA/2-1)/NfftA * PRF;

D                               = sqrt(1-(lambda*Fnan/2/Vr).^2); % 徙动因子[式(6.24)]
Dc                              = sqrt(1-(lambda*Fdc/(2*Vr))^2);
% WinA                            = kaiser(NfftA, 2.5); % 生成Kaiser窗
WinA                            = ones(NfftA, 1);

Ha                              = exp(1j*4*pi*Rrefs.'*D*fc/C) .* WinA.'; % 方位向匹配滤波器[式(6.26)]
eta_c                           = -lambda*Rrefs.'*Fdc/(2*Vr^2*Dc);
Ha                              = Ha .* exp(-1j*2*pi*eta_c*fbase);
temp                            = fftshift(ifft(ifftshift(DATA,2),[],2),2); % RD域变换回时域
temp                            = padarray(temp, [0 (NfftA - Na)/2], 0); % 方位向两边补零
clear DATA;

% temp = fftshift(ifft(...                            % 方位向匹配滤波
%                          fft(ifftshift(temp, 2),[],2) .* ifftshift(Ha, 2),[],2), 2);

for rr = 1:Nr
    if mod(rr, 128) == 0
        fprintf('AC Processing: %d / %d\n', rr, Nr);
    end
    temp(rr,:) = fftshift(ifft(fft(ifftshift(temp(rr,:))).*ifftshift(Ha(rr,:))));
end

DATA                            = temp(:,NfftA / 2 - Na / 2 + 1 : NfftA / 2 + Na / 2);
clear D Ha temp;

%% 显示图像
figure;imagesc(abs(DATA));

img                             = abs(DATA/max(DATA,[],'all'));
figure;
imagesc(img,[0,0.02]);
colormap('gray');

figure;
imagesc(db(img));
colormap('gray');

Imglog      = mat2gray(log2(img));
ImgGamma    = imadjust(Imglog,[],[],4);
figure;imshow(flip(Imglog));
figure;imshow(flip(ImgGamma));
imwrite(single(flip(ImgGamma)),'RD_SM.jpg');