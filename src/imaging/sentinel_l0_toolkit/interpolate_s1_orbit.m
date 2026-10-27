function info = interpolate_s1_orbit(info)
%INTERPOLATE_S1_ORBIT Place ECEF position/velocity on every selected PRI row.
%   info = interpolate_s1_orbit(info) preserves the original navigation records
%   and adds pulse_position_ecef_m / pulse_velocity_ecef_mps to info.orbit.
%   Target times are info.timeline.time_relative_s, measured from the FIRST
%   SOURCE PACKET, including calibration/noise/missing rows. StartPulse is not
%   added again. Native orbit times must use the same source origin.
%   A piecewise cubic Hermite curve uses native positions and velocities;
%   interpolated velocity is the derivative of that same curve. Outside the
%   native time range, constant endpoint velocity is allowed for at most one
%   adjacent native navigation interval, with an explicit extrapolation mask.
%   All absolute times are GPS seconds, not UTC. Calculations use relative
%   seconds to avoid subtracting large absolute times during interpolation.

assert(isstruct(info) && isscalar(info), 'S1Orbit:InvalidInput', 'info must be a scalar struct.');
assert(isfield(info, 'timeline') && isstruct(info.timeline) && ...
    isfield(info.timeline, 'time_relative_s'), 'S1Orbit:MissingTimeline', ...
    'info.timeline.time_relative_s is required.');
target = info.timeline.time_relative_s;
validateattributes(target, {'numeric'}, {'vector','real','finite'}, mfilename, 'timeline.time_relative_s');
target = double(target(:));
assert(all(diff(target) > 0), 'S1Orbit:InvalidTimeline', ...
    'Selected PRI times must be strictly increasing.');
if ~isfield(info, 'orbit') || isempty(info.orbit)
    info.orbit = struct('available', false);
end
assert(isstruct(info.orbit) && isscalar(info.orbit), 'S1Orbit:InvalidOrbit', ...
    'info.orbit must be a scalar struct.');
orbit = info.orbit;
n = numel(target);
orbit.pulse_time_relative_s = target;
orbit.pulse_time_gps_s = NaN(n, 1);
orbit.pulse_position_ecef_m = NaN(n, 3);
orbit.pulse_velocity_ecef_mps = NaN(n, 3);
orbit.pulse_extrapolated_mask = false(n, 1);
orbit.pulse_extrapolated_count = 0;
orbit.pulse_count = n;
orbit.pulse_interpolation_available = false;
orbit.interpolation_method = 'Piecewise cubic Hermite constrained by native ECEF position and velocity; velocity is its derivative';
orbit.extrapolation_note = ['Constant endpoint velocity outside the native interval; limited to one adjacent ' ...
    'native navigation interval at each end. Extrapolated rows are marked, including zero PRI slots.'];
orbit.pulse_time_note = ['GPS seconds, NOT UTC. Relative times use the first source packet as origin; ' ...
    'the selected PRI offset is already present and is not added again.'];
orbit.pulse_units_note = 'Position: ECEF metres; velocity: ECEF metres/second; time: seconds';
orbit.pulse_interpolation_reason = '';
orbit.extrapolation_limit_before_s = NaN;
orbit.extrapolation_limit_after_s = NaN;

if isfield(orbit, 'available')
    assert((islogical(orbit.available) || isnumeric(orbit.available)) && ...
        isscalar(orbit.available) && isfinite(double(orbit.available)) && ...
        ismember(double(orbit.available), [0, 1]), 'S1Orbit:InvalidOrbit', ...
        'orbit.available must be a scalar logical.');
    if ~orbit.available
        info.orbit = unavailable(orbit, 'Native orbit is unavailable; per-PRI ECEF arrays remain NaN.');
        return
    end
end
required = {'time_relative_s','time_gps_s','time_origin_gps_s','position_ecef_m','velocity_ecef_mps'};
if ~all(isfield(orbit, required))
    info.orbit = unavailable(orbit, 'Complete native orbit time/position/velocity fields are unavailable.');
    return
end
validateattributes(orbit.time_origin_gps_s, {'numeric'}, {'scalar','real'}, mfilename, 'orbit.time_origin_gps_s');
origin = double(orbit.time_origin_gps_s);
if isnan(origin)
    info.orbit = unavailable(orbit, 'The shared first-packet GPS origin is unknown; per-PRI ECEF alignment is unavailable.');
    return
end
assert(isfinite(origin), 'S1Orbit:InvalidTime', 'The GPS time origin must be finite.');
orbit.pulse_time_gps_s = origin + target;
validateattributes(orbit.time_relative_s, {'numeric'}, {'vector','real','finite'}, mfilename, 'orbit.time_relative_s');
validateattributes(orbit.time_gps_s, {'numeric'}, {'vector','real','finite'}, mfilename, 'orbit.time_gps_s');
nativeTime = double(orbit.time_relative_s(:));
nativeGps = double(orbit.time_gps_s(:));
m = numel(nativeTime);
assert(numel(nativeGps) == m && all(diff(nativeTime) > 0), 'S1Orbit:InvalidTime', ...
    'Native navigation times must be strictly increasing with matching GPS/relative lengths.');
validateattributes(orbit.position_ecef_m, {'numeric'}, {'2d','real','finite'}, mfilename, 'orbit.position_ecef_m');
validateattributes(orbit.velocity_ecef_mps, {'numeric'}, {'2d','real','finite'}, mfilename, 'orbit.velocity_ecef_mps');
assert(isequal(size(orbit.position_ecef_m), [m, 3]) && ...
    isequal(size(orbit.velocity_ecef_mps), [m, 3]), 'S1Orbit:InvalidShape', ...
    'Native position and velocity must both be N x 3, with N matching the time vector.');
if m < 2
    info.orbit = unavailable(orbit, 'At least two native position/velocity records are required.');
    return
end
if isfield(orbit, 'reference_frame')
    assert(strcmpi(strtrim(char(orbit.reference_frame)), 'ECEF'), 'S1Orbit:InvalidFrame', ...
        'Native position and velocity must use ECEF coordinates.');
end
if isfield(orbit, 'time_system')
    assert(strcmpi(strtrim(char(orbit.time_system)), 'GPS'), 'S1Orbit:InvalidTimeSystem', ...
        'Native absolute times must use GPS seconds, not UTC.');
end
% Absolute GPS additions cannot resolve much below eps(~1e9) seconds. This
% tolerance is used ONLY to validate provenance, never to compute the curve.
absoluteTolerance = 8 * eps(max([1; abs(origin); abs(nativeGps)]));
assert(all(abs(nativeGps - (origin + nativeTime)) <= absoluteTolerance), ...
    'S1Orbit:TimeOriginMismatch', 'Native GPS/relative times disagree with orbit.time_origin_gps_s.');
if isfield(info.timeline, 'time_origin_gps_s')
    checkOrigin(info.timeline.time_origin_gps_s, origin, absoluteTolerance, 'timeline');
end
if isfield(info, 'packet_metadata') && isstruct(info.packet_metadata) && ...
        isfield(info.packet_metadata, 'values')
    packet = info.packet_metadata.values;
    if isstruct(packet) && all(isfield(packet, {'Coarse_Time','Fine_Time'})) && ...
            ~isempty(packet.Coarse_Time) && ~isempty(packet.Fine_Time)
        checkOrigin(double(packet.Coarse_Time(1)) + double(packet.Fine_Time(1)), ...
            origin, absoluteTolerance, 'first source packet');
    end
end
orbit.time_alignment_note = ['Native GPS/relative times checked against their origin; available timeline/first-packet ' ...
    'origins checked against it. All interpolation uses the supplied relative times.'];

position = double(orbit.position_ecef_m);
velocity = double(orbit.velocity_ecef_mps);
intervals = diff(nativeTime);
before = target < nativeTime(1);
after = target > nativeTime(end);
inside = ~(before | after);
relativeTolerance = 32 * eps(max([1; abs(nativeTime); abs(target)]));
assert(~any(nativeTime(1) - target(before) > intervals(1) + relativeTolerance) && ...
    ~any(target(after) - nativeTime(end) > intervals(end) + relativeTolerance), ...
    'S1Orbit:ExtrapolationTooLong', ...
    'Requested PRI times exceed the permitted one-native-interval endpoint extrapolation.');
orbit.extrapolation_limit_before_s = intervals(1);
orbit.extrapolation_limit_after_s = intervals(end);

if any(inside)
    query = target(inside);
    segment = discretize(query, nativeTime);
    h = intervals(segment);
    dt = query - nativeTime(segment);
    p0 = position(segment, :);
    p1 = position(segment + 1, :);
    v0 = velocity(segment, :);
    v1 = velocity(segment + 1, :);
    delta = p1 - p0;
    c2 = (3 * delta ./ h - 2 * v0 - v1) ./ h;
    c3 = (-2 * delta ./ h + v0 + v1) ./ (h .^ 2);
    orbit.pulse_position_ecef_m(inside, :) = p0 + dt .* (v0 + dt .* (c2 + dt .* c3));
    orbit.pulse_velocity_ecef_mps(inside, :) = v0 + dt .* (2 * c2 + 3 * dt .* c3);
end
if any(before)
    orbit.pulse_position_ecef_m(before, :) = position(1, :) + (target(before) - nativeTime(1)) .* velocity(1, :);
    orbit.pulse_velocity_ecef_mps(before, :) = repmat(velocity(1, :), sum(before), 1);
end
if any(after)
    orbit.pulse_position_ecef_m(after, :) = position(end, :) + (target(after) - nativeTime(end)) .* velocity(end, :);
    orbit.pulse_velocity_ecef_mps(after, :) = repmat(velocity(end, :), sum(after), 1);
end
% Preserve stored state values exactly at knots, including the final endpoint.
[atKnot, nativeRow] = ismember(target, nativeTime);
orbit.pulse_position_ecef_m(atKnot, :) = position(nativeRow(atKnot), :);
orbit.pulse_velocity_ecef_mps(atKnot, :) = velocity(nativeRow(atKnot), :);
assert(all(isfinite(orbit.pulse_position_ecef_m), 'all') && ...
    all(isfinite(orbit.pulse_velocity_ecef_mps), 'all'), 'S1Orbit:InvalidResult', ...
    'Orbit interpolation produced non-finite values.');
orbit.pulse_extrapolated_mask = before | after;
orbit.pulse_extrapolated_count = sum(orbit.pulse_extrapolated_mask);
orbit.pulse_interpolation_available = true;
info.orbit = orbit;
end

function orbit = unavailable(orbit, reason)
orbit.available = false;
orbit.pulse_interpolation_reason = reason;
end

function checkOrigin(candidate, expected, tolerance, source)
validateattributes(candidate, {'numeric'}, {'scalar','real','finite'});
assert(abs(double(candidate) - expected) <= tolerance, 'S1Orbit:TimeOriginMismatch', ...
    'GPS origin from %s differs from the native orbit origin.', source);
end
