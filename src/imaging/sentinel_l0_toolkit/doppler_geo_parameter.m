function info = doppler_geo_parameter(info)
%DOPPLER_GEO_PARAMETER 在 Shift_Delay_correction 后读取多普勒和计算几何参数。
% info = doppler_geo_parameter(info)，不改变 echo。
% info.doppler：当前选段的一组距离向多普勒中心、方位调频率及其 XML 来源。
% info.xml_orbit：XML 原始轨道及逐 PRI 插值轨道，不覆盖 info.orbit。
% info.geometry：参考距离、目标位置、卫星速率、逐距离门等效速度等。
% 参考 MATLAB：parse_xmlpara/get_doppler_parameter/get_target_pos/get_equivalent_vel。
% 定位沿用源码的零多普勒、WGS84 零高程、右视模型，不是 DEM 精确定位。

assert(isfield(info, 'shift_delay_correction') && info.shift_delay_correction.applied, ...
    'doppler_geo_parameter:ProcessingOrder', '请先运行 Shift_Delay_correction。');
correction = info.shift_delay_correction;
xml = readstruct(correction.xml_file);                  % 复用已经匹配好的 SLC XML
time = double(info.timeline.time_relative_s(:));       % 完整 PRI 时间轴，包含零行
pulseUtc = double(correction.pulse_time_posix_s(:));
rangeTime = double(info.radar.slant_range_time_s(:));   % 双程时间，含 Rank*PRI
range = double(info.radar.slant_range_m(:));
Na = numel(time); Nr = numel(range);
assert(Na >= 3 && numel(pulseUtc) == Na && all(isfinite([time;pulseUtc])) && ...
    all(diff(time) > 0), 'doppler_geo_parameter:InvalidTime', '等效速度计算至少需要 3 个有序脉冲时刻。');
assert(Nr >= 2 && numel(rangeTime) == Nr && all(isfinite([range;rangeTime])) && ...
    all(range > 0) && all(diff(range) > 0), ...
    'doppler_geo_parameter:InvalidRange', '请使用 echo_SWST 返回的完整距离轴。');
centerA = floor(Na/2)+1;
centerR = double(info.radar.center_sample_index);       % 偶数 Nr 时为 Nr/2+1
assert(centerR == floor(Nr/2)+1, 'doppler_geo_parameter:CenterIndex', '中心距离采样点与 Nr 不一致。');
lambda = double(info.radar.lambda);
prf = double(info.radar.prf_hz);
validateattributes(lambda, {'double'}, {'scalar','finite','positive'});
validateattributes(prf, {'double'}, {'scalar','finite','positive'});
processing = xml.imageAnnotation.processingInformation;
product = xml.generalAnnotation.productInformation;

%% 1. XML 轨道：UTC、ECEF、m、m/s，按源码做 spline 插值
orbit = xml.generalAnnotation.orbitList.orbit;
orbitTime = utcSeconds([orbit.time].');
assert(numel(orbitTime) >= 2 && all(diff(orbitTime) > 0) && ...
    all(strcmpi(string({orbit.frame}), 'Earth Fixed')), ...
    'doppler_geo_parameter:XMLOrbit', 'XML 轨道应为有序的 Earth Fixed 状态向量。');
position = xyz([orbit.position]);
velocity = xyz([orbit.velocity]);
assert(all(isfinite([position,velocity]), 'all'), ...
    'doppler_geo_parameter:XMLOrbit', 'XML 轨道位置与速度必须为有限值。');
% 插值使用相对秒数，保持 PRI 间隔精度；绝对时刻来自上一步 GPS→UTC 转换。
queryTime = pulseUtc(1)-orbitTime(1) + time-time(1);
xmlPosition = interp1(orbitTime-orbitTime(1), position, queryTime, 'spline', 'extrap');
xmlVelocity = interp1(orbitTime-orbitTime(1), velocity, queryTime, 'spline', 'extrap');
info.xml_orbit = struct('xml_file', correction.xml_file, 'reference_frame', 'ECEF', ...
    'time_posix_s', orbitTime, 'position_ecef_m', position, 'velocity_ecef_mps', velocity, ...
    'pulse_time_posix_s', pulseUtc, 'pulse_position_ecef_m', xmlPosition, ...
    'pulse_velocity_ecef_mps', xmlVelocity, ...
    'extrapolated_mask', queryTime < 0 | queryTime > orbitTime(end)-orbitTime(1), ...
    'method', 'Independent position/velocity spline interpolation, as in parse_xmlpara');

% 与源码一样优先使用 RAW 轨道定位；RAW 逐脉冲轨道不可用时明确记录 XML 回退。
satPosition = xmlPosition; satVelocity = xmlVelocity;
orbitSource = 'XML';
if isfield(info, 'orbit') && all(isfield(info.orbit, ...
        {'pulse_position_ecef_m','pulse_velocity_ecef_mps'}))
    p = double(info.orbit.pulse_position_ecef_m);
    v = double(info.orbit.pulse_velocity_ecef_mps);
    if isequal(size(p), [Na,3]) && isequal(size(v), [Na,3]) && all(isfinite([p,v]), 'all')
        satPosition = p; satVelocity = v; orbitSource = 'RAW';
    end
end
satelliteSpeed = vecnorm(satVelocity, 2, 2);
assert(all(isfinite(satPosition), 'all') && all(isfinite(satelliteSpeed) & satelliteSpeed > 0), ...
    'doppler_geo_parameter:InvalidOrbit', '用于几何计算的轨道无效。');

%% 2. 多普勒：与源码一样，为当前选段各选一条 DC/FM 记录，再沿距离求值
dcList = xml.dopplerCentroid.dcEstimateList.dcEstimate;
fmList = xml.generalAnnotation.azimuthFmRateList.azimuthFmRate;
[~, dcIndex] = min(sum(abs(pulseUtc - utcSeconds([dcList.azimuthTime].').'), 1));
[~, fmIndex] = min(sum(abs(pulseUtc - utcSeconds([fmList.azimuthTime].').'), 1));
dc = dcList(dcIndex); fm = fmList(fmIndex);
geometryCoefficients = polynomial(dc.geometryDcPolynomial);
dataCoefficients = polynomial(dc.dataDcPolynomial);
if isfield(fm, 'azimuthFmRatePolynomial')
    fmCoefficients = polynomial(fm.azimuthFmRatePolynomial);
else
    fmCoefficients = double([fm.c0, fm.c1, fm.c2]);      % 兼容早期 XML
end
geometryDc = polyval(fliplr(geometryCoefficients), rangeTime-double(dc.t0));
dataDc = polyval(fliplr(dataCoefficients), rangeTime-double(dc.t0));
fmRate = polyval(fliplr(fmCoefficients), rangeTime-double(fm.t0));
assert(all(isfinite([geometryDc;dataDc;fmRate])), ...
    'doppler_geo_parameter:Polynomial', 'XML 多普勒多项式求值结果无效。');
referenceDc = fix(geometryDc/prf)*prf + dataDc;         % 保留原源码公式供对照
% XML dataDcPolynomial 已是绝对多普勒，不能再次加 PRF 整周数。
badDataDc = strcmpi(string(dc.dataDcRmsErrorAboveThreshold), 'true');
centroid = dataDc; centroidSource = 'XML dataDcPolynomial (absolute Doppler)';
if badDataDc
    centroid = geometryDc; centroidSource = 'XML geometryDcPolynomial (data DC failed quality flag)';
end
steeringRate = deg2rad(double(product.azimuthSteeringRate));
swathProcessing = processing.swathProcParamsList.swathProcParams;
swathProcessing = swathProcessing(strcmpi(string({swathProcessing.swath}), string(xml.adsHeader.swath)));
assert(isscalar(swathProcessing), 'doppler_geo_parameter:Swath', 'XML 必须包含当前子带唯一的处理参数。');
info.doppler = struct('xml_file', correction.xml_file, ...
    'centroid_hz', centroid, 'geometry_centroid_hz', geometryDc, 'data_centroid_hz', dataDc, ...
    'reference_formula_centroid_hz', referenceDc, 'centroid_source', centroidSource, ...
    'azimuth_fm_rate_hz_s', fmRate, ...
    'centroid_rate_hz_s', 2*satelliteSpeed(centerA)/lambda*steeringRate, ...
    'processing_bandwidth_hz', double(swathProcessing.azimuthProcessing.processingBandwidth), ...
    'dc_estimate_index', dcIndex, 'fm_estimate_index', fmIndex, ...
    'dc_estimate_time_posix_s', utcSeconds(dc.azimuthTime), ...
    'fm_estimate_time_posix_s', utcSeconds(fm.azimuthTime), ...
    'dc_t0_s', double(dc.t0), 'fm_t0_s', double(fm.t0), ...
    'geometry_coefficients_ascending', geometryCoefficients, ...
    'data_coefficients_ascending', dataCoefficients, 'fm_coefficients_ascending', fmCoefficients, ...
    'data_rms_error_hz', double(dc.dataDcRmsError), 'data_rms_error_above_threshold', badDataDc, ...
    'center_centroid_hz', centroid(centerR), 'center_azimuth_fm_rate_hz_s', fmRate(centerR), ...
    'xml_dc_estimates', dcList, 'xml_fm_estimates', fmList, ...
    'note', 'One representative DC/FM record for this selected segment; bandwidth is XML processing bandwidth.');

%% 3. 每隔 30 个距离门定位，拟合 R(t)^2 = R0^2 + Veq^2*(t-t0)^2
a = 6378137; b = a*(1-1/298.257223563);               % WGS84，与参考源码一致
sampleIndices = unique([1:30:Nr, centerR, Nr]).';      % 加入中心和末列，不在末端外推
count = numel(sampleIndices);
targets = zeros(count,3); residuals = zeros(count,3);
equivalentSpeed = zeros(count,1); closestRange = zeros(count,1); closestRow = zeros(count,1);
for k = 1:count
    [target, residuals(k,:)] = locateTarget(satPosition(centerA,:), ...
        satVelocity(centerA,:), range(sampleIndices(k)), a, b);
    targets(k,:) = target;
    distanceSquared = sum((satPosition-target).^2, 2);
    [minimumSquared, row] = min(distanceSquared);
    eta = time-time(row);
    displacement = satPosition-satPosition(row,:);
    % 等价于 R(t)^2-R0^2，避免两个约 1e12 的大数直接相减。
    deltaRangeSquared = sum(displacement.^2,2) + 2*displacement*(satPosition(row,:)-target).';
    speedSquared = sum(deltaRangeSquared.*eta.^2)/sum(eta.^4);
    assert(isfinite(speedSquared) && speedSquared > 0, ...
        'doppler_geo_parameter:EquivalentVelocity', '等效速度拟合无效，请检查轨道及所选时间区间。');
    equivalentSpeed(k) = sqrt(speedSquared);
    closestRange(k) = sqrt(minimumSquared); closestRow(k) = row;
end
velocityByRange = interp1(sampleIndices, equivalentSpeed, (1:Nr).', 'spline');
assert(all(isfinite(velocityByRange) & velocityByRange > 0), ...
    'doppler_geo_parameter:EquivalentVelocity', '等效速度插值无效。');
centerTarget = targets(sampleIndices == centerR,:);
latitude = atan2d(centerTarget(3)*a^2, hypot(centerTarget(1),centerTarget(2))*b^2);
longitude = atan2d(centerTarget(2), centerTarget(1));

%% 4. 存储几何参数；XML 参考距离与当前中心距离分别保存
info.geometry = struct('xml_file', correction.xml_file, 'orbit_source', orbitSource, ...
    'xml_reference_range_m', double(processing.referenceRange), ...
    'xml_first_slant_range_time_s', double(xml.imageAnnotation.imageInformation.slantRangeTime), ...
    'xml_radar_frequency_hz', double(product.radarFrequency), ...
    'azimuth_steering_rate_rad_s', steeringRate, 'satellite_speed_mps', satelliteSpeed, ...
    'center_azimuth_index', centerA, 'center_range_index', centerR, ...
    'center_time_posix_s', pulseUtc(centerA), 'center_slant_range_m', range(centerR), ...
    'target_position_ecef_m', centerTarget, 'target_lla_deg_m', [latitude,longitude,0], ...
    'sample_range_indices', sampleIndices, 'target_sample_position_ecef_m', targets, ...
    'target_residuals_m', residuals, 'equivalent_velocity_mps', velocityByRange, ...
    'equivalent_velocity_sample_mps', equivalentSpeed, ...
    'center_equivalent_velocity_mps', velocityByRange(centerR), ...
    'closest_range_sample_m', closestRange, 'closest_time_sample_posix_s', pulseUtc(closestRow), ...
    'look_direction', -1, 'target_height_m', 0, 'ellipsoid_axes_m', [a,b], ...
    'model_note', 'Zero-Doppler, zero-height, right-looking WGS84 model; not a DEM/geometric calibration.', ...
    'reference_range_note', 'XML referenceRange is for radiometric range-spreading-loss normalization, not scene-center range.');
if isfield(xml.imageAnnotation.imageInformation, 'zeroDopMinusAcqTime')
    info.geometry.xml_zero_doppler_minus_acquisition_s = double(xml.imageAnnotation.imageInformation.zeroDopMinusAcqTime);
end
% 源码固定配置值，明确与 XML 实测/读取值分开。
info.geometry.reference_defaults = struct('antenna_azimuth_length_m', 12.3, ...
    'antenna_elevation_length_m', 0.82, 'beam_azimuth_width_rad', deg2rad(0.23), ...
    'beam_elevation_width_rad', deg2rad(3.43), 'azimuth_steering_extent_rad', deg2rad([-0.9,0.9]), ...
    'elevation_steering_extent_rad', deg2rad([-13,12.3]));
% 波束地址换算的是指令角，不是由多普勒中心反算的斜视角。
if isfield(info, 'packet_metadata') && isfield(info.packet_metadata, 'values')
    p = info.packet_metadata.values;
    if all(isfield(p, {'Azimuth_Beam_Address','Signal_Type','selected_row_zero_based'}))
        selected = p.Signal_Type(:) == 0 & p.selected_row_zero_based(:) >= 0 & p.selected_row_zero_based(:) < Na;
        angle = nan(Na,1);
        angle(double(p.selected_row_zero_based(selected))+1) = ...
            deg2rad(double(p.Azimuth_Beam_Address(selected))/1023*1.8-0.9);
        info.geometry.beam_command_squint_rad = angle;
    end
end
end

function [target, residual] = locateTarget(satellite, velocity, range, a, b)
% 解三个方程：斜距球面、零多普勒平面、WGS84 椭球。
s = satellite(:); v = velocity(:)/norm(velocity);
axesSquared = [a;a;b].^2;
earthRadius = 1/sqrt(sum((s/norm(s)).^2./axesSquared));
circleCenter = dot(s,v)*v;
direction = circleCenter-s; distance = norm(direction);
x = (range^2-earthRadius^2+dot(s,v)^2+distance^2)/(2*distance);
assert(distance > 0 && range^2 > x^2, 'doppler_geo_parameter:NoIntersection', ...
    '此斜距没有有效的右视地表交点，请核对距离轴。');
direction = direction/distance;
target = s + x*direction - sqrt(range^2-x^2)*cross(v,direction); % 选择右视解
for iteration = 1:20
    delta = target-s; r = norm(delta);
    residual = [r-range; dot(v,delta); (sum(target.^2./axesSquared)-1)*a/2];
    if max(abs(residual)) < 1e-4, break; end
    jacobian = [delta.'/r; v.'; (a*target./axesSquared).'];
    target = target-jacobian\residual;
end
delta = target-s;
residual = [norm(delta)-range; dot(v,delta); (sum(target.^2./axesSquared)-1)*a/2];
assert(all(isfinite(target)) && max(abs(residual)) < 1e-3 && dot(delta,cross(v,s)) > 0, ...
    'doppler_geo_parameter:GeolocationFailed', '目标定位未收敛到右视解。');
target = target.'; residual = residual.';
end

function values = polynomial(field)
% XML 按常数项、一次项、二次项排列；polyval 调用时再翻转。
if isstruct(field), field = field.Text; end
if isnumeric(field), values = double(field(:).');
else, values = sscanf(char(field), '%f').'; end
assert(~isempty(values) && all(isfinite(values)), ...
    'doppler_geo_parameter:Polynomial', 'XML 多项式系数无效。');
end

function time = utcSeconds(value)
if ~isdatetime(value)
    value = datetime(value, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");
end
value.TimeZone = 'UTC';
time = posixtime(value(:));
assert(all(isfinite(time)), 'doppler_geo_parameter:XMLTime', 'XML 时刻无效。');
end

function values = xyz(records)
values = double([[records.x].', [records.y].', [records.z].']);
end
