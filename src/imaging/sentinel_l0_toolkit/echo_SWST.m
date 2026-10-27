function [echoAligned, info] = echo_SWST(echo, info)
%ECHO_SWST Align RAW receive windows to the minimum ECHO SWST.
%   [echoAligned, infoAligned] = echo_SWST(echo, info)
%   Rows are the complete PRI timeline returned by read_s1_raw.
%   Positive SWST differences delay data to the right by a frequency-domain
%   phase ramp. Append ceil(max SWST difference * Fs) zeros, then make Nr even.
%   Numeric single/double input returns a numeric complex matrix of that class.
%   A matfile input exposing iq returns a NEW matfile, processed in row blocks.
%   timeline.swst_s becomes the common origin; the original per-row values
%   are retained in swst_alignment.original_swst_s.
%   radar.slant_range_m is an uncalibrated RAW timing axis including Rank*PRI
%   and the Sentinel-1 decimation transient suppression delay, not a focused
%   or geometrically calibrated range axis. The centre index is Nr/2+1.

diskInput = isa(echo, 'matlab.io.MatFile');
if diskInput
    variables = whos(echo);
    item = variables(strcmp({variables.name}, 'iq'));
    assert(numel(item) == 1, 'echo_SWST:MissingIQ', 'The matfile must contain iq.');
    shape = size(echo, 'iq');
    precision = item.class;
else
    assert(isnumeric(echo) && ismatrix(echo) && ~issparse(echo), ...
        'echo_SWST:InvalidEcho', 'echo must be a dense numeric matrix or a matfile containing iq.');
    shape = size(echo);
    precision = class(echo);
end
assert(numel(shape) == 2 && all(shape > 0) && ismember(precision, {'single','double'}), ...
    'echo_SWST:InvalidEcho', 'Expected a nonempty two-dimensional single/double matrix.');
Na = shape(1);
NrInput = shape(2);
assert(isstruct(info) && isscalar(info) && isfield(info, 'timeline') && isfield(info, 'radar'), ...
    'echo_SWST:MissingInfo', 'info.timeline and info.radar are required.');
assert(isfield(info.timeline, 'echo_mask') && isfield(info.timeline, 'swst_s'), ...
    'echo_SWST:MissingInfo', 'The per-row echo_mask and swst_s are required.');
echoMask = info.timeline.echo_mask(:);
swst = double(info.timeline.swst_s(:));
assert(numel(echoMask) == Na && numel(swst) == Na && ...
    all(ismember(echoMask, [0,1])), 'echo_SWST:RowMismatch', ...
    'echo rows and info.timeline rows must agree; subset the metadata together with the data.');
echoMask = logical(echoMask);
if isfield(info, 'swst_alignment') && isfield(info.swst_alignment, 'applied') && info.swst_alignment.applied
    assert(isfield(info.radar, 'Nr') && info.radar.Nr == NrInput && ...
        numel(info.swst_alignment.shift_samples) == Na, 'echo_SWST:AlreadyAlignedShape', ...
        'The already-aligned metadata does not match this matrix.');
    echoAligned = echo;
    return
end
assert(any(echoMask), 'echo_SWST:NoEchoRows', 'No ECHO rows are available to define the SWST origin.');
assert(all(isfinite(swst(echoMask)) & swst(echoMask) >= 0), ...
    'echo_SWST:InvalidSWST', 'Every ECHO row needs a finite, nonnegative SWST.');
assert(isfield(info.radar, 'range_sample_rate_hz') && isfield(info.radar, 'pri_s'), ...
    'echo_SWST:MissingInfo', 'range_sample_rate_hz and pri_s are required.');
Fs = double(info.radar.range_sample_rate_hz);
PRI = double(info.radar.pri_s);
assert(isscalar(Fs) && isfinite(Fs) && Fs > 0 && isscalar(PRI) && isfinite(PRI) && PRI > 0, ...
    'echo_SWST:InvalidTiming', 'Fs and PRI must be finite positive scalars.');
rank = getRank(info, Na);
swstMin = min(swst(echoMask));
swstMax = max(swst(echoMask));
shiftSamples = zeros(Na, 1);
shiftSamples(echoMask) = (swst(echoMask)-swstMin)*Fs;
% Snap only numerical roundoff near integer shifts, before using ceil.
integerTolerance = 32*eps(max(1, max(abs(swst(echoMask)))*Fs));
nearest = round(shiftSamples);
nearInteger = abs(shiftSamples-nearest) <= integerTolerance;
shiftSamples(nearInteger) = nearest(nearInteger);
shiftPadding = ceil(max(shiftSamples));
Nr = NrInput + shiftPadding;
evenPadding = mod(Nr, 2);
Nr = Nr + evenPadding;
assert(isfinite(Nr) && Nr <= flintmax, 'echo_SWST:InvalidLength', 'The padded range length is invalid.');

info.swst_alignment = struct('applied', true, 'original_swst_s', swst, ...
    'row_applied_mask', echoMask, 'reference_swst_s', swstMin, 'maximum_swst_s', swstMax, ...
    'shift_samples', shiftSamples, 'shift_s', shiftSamples/Fs, ...
    'original_range_samples', NrInput, 'padding_samples', Nr-NrInput, ...
    'shift_padding_samples', shiftPadding, 'even_padding_samples', evenPadding, ...
    'integer_snap_tolerance_samples', integerTolerance, ...
    'method', 'FFT along range; exp(-1i*2*pi*f*deltaSWST); inverse FFT', ...
    'boundary_note', 'Input is zero-padded before shifting. Fractional shifts use periodic bandlimited interpolation and may have sinc tails.');
if isfield(info, 'mat_file'), info.swst_alignment.source_mat_file = info.mat_file; end
if isfield(info.radar, 'swl_s'), info.swst_alignment.original_swl_s = info.radar.swl_s; end
if isfield(info.timeline, 'range_window_offset_samples')
    info.swst_alignment.original_range_window_offset_samples = info.timeline.range_window_offset_samples;
end
sampleTime = (0:Nr-1)'/Fs;
referenceFrequency = 37.53472224e6;
transientDelay = 320/(8*referenceFrequency);
c = 299792458;
rangeTime = rank*PRI + swstMin + transientDelay + sampleTime;
info.timeline.swst_s = repmat(swstMin, Na, 1);
info.timeline.aligned_swst_s = info.timeline.swst_s;
info.timeline.range_window_offset_samples = zeros(Na, 1);
info.timeline.range_window_offset_note = 'All rows use the common minimum ECHO SWST after alignment; non-ECHO rows remain zero.';
info.radar.Nr = Nr;
info.radar.swst_s = swstMin;
info.radar.swst_reference_s = swstMin;
info.radar.shared_absolute_range_origin = true;
info.radar.rank = rank;
info.radar.range_sample_interval_s = 1/Fs;
info.radar.range_sample_time_relative_s = sampleTime;
info.radar.sampling_window_time_s = swstMin + sampleTime;
info.radar.slant_range_time_s = rangeTime;
info.radar.slant_range_m = c/2*rangeTime;
info.radar.range_sample_spacing_m = c/(2*Fs);
info.radar.aligned_window_duration_s = Nr/Fs;
info.radar.center_sample_index = Nr/2+1;
info.radar.center_slant_range_m = info.radar.slant_range_m(info.radar.center_sample_index);
info.radar.reference_frequency_hz = referenceFrequency;
info.radar.decimation_transient_delay_s = transientDelay;
info.radar.swst_bias_applied = false;
info.radar.swst_bias_s = NaN;
info.radar.range_axis_note = ['Common uncalibrated RAW timing axis: c/2*(Rank*PRI + minSWST + ' ...
    '320/(8*f_ref) + sample/Fs). No pulse compression, matched-filter reference shift, or SWST-bias correction.'];
if isfield(info.radar, 'parameter_sources')
    info.radar.parameter_sources.swst_s = 'Minimum selected ECHO SWST after alignment; original values in swst_alignment.original_swst_s';
end
if isfield(info, 'processing_info')
    info.processing_info.range_samples = Nr;
    info.processing_info.processing_level = 'Level-0 decoded I/Q with SWST alignment; no range/azimuth compression';
end

% Never alter an input MAT file. A disk-backed result is published only after
% every row and the complete metadata have been written to a new file.
if diskInput
    sourcePath = echo.Properties.Source;
    [outputDir, stem] = fileparts(sourcePath);
    if isempty(outputDir), outputDir = pwd; end
    [~, token] = fileparts(tempname(outputDir));
    finalPath = fullfile(outputDir, [stem '_SWST_' token '.mat']);
    temporaryPath = [finalPath '.partial.mat'];
    cleanup = onCleanup(@()removeTemporary(temporaryPath));
    iq = complex(zeros(1, 1, precision)); %#ok<NASGU>
    save(temporaryPath, 'iq', '-v7.3');
    writer = matfile(temporaryPath, 'Writable', true);
    writer.iq(Na, Nr) = complex(cast(0, precision));
    info.mat_file = finalPath;
    info.swst_alignment.output_file = finalPath;
else
    echoAligned = complex(zeros(Na, Nr, precision));
    info.mat_file = '';
    info.swst_alignment.output_file = '';
end

bins = [0:Nr/2-1, -Nr/2:-1];
batchRows = 64;
for first = 1:batchRows:Na
    rows = first:min(first+batchRows-1, Na);
    if diskInput, sourceBlock = echo.iq(rows, :); else, sourceBlock = echo(rows, :); end
    assert(all(isfinite(sourceBlock), 'all'), 'echo_SWST:NonfiniteEcho', 'Echo data must be finite.');
    active = echoMask(rows);
    assert(all(sourceBlock(~active, :) == 0, 'all'), 'echo_SWST:NonzeroInactiveRow', ...
        'Non-ECHO rows must be zero; the mask must describe this echo matrix.');
    block = complex(zeros(numel(rows), Nr, precision));
    localShifts = shiftSamples(rows);
    unchanged = active & localShifts == 0;
    block(unchanged, 1:NrInput) = sourceBlock(unchanged, :);
    for delay = unique(localShifts(active & localShifts ~= 0))'
        selected = active & localShifts == delay;
        phase = cast(exp(-2i*pi*(delay/Nr)*bins), precision);
        block(selected, :) = ifft(fft(sourceBlock(selected, :), Nr, 2).*phase, [], 2);
    end
    if diskInput, writer.iq(rows, :) = block; else, echoAligned(rows, :) = block; end
end
if diskInput
    clear writer
    save(temporaryPath, '-struct', 'info', '-append');
    assert(~isfile(finalPath), 'echo_SWST:OutputExists', 'The new output path unexpectedly exists.');
    [ok, message] = movefile(temporaryPath, finalPath);
    assert(ok, 'echo_SWST:PublishFailed', '%s', message);
    clear cleanup
    echoAligned = matfile(finalPath);
end
end

function rank = getRank(info, Na)
if isfield(info.radar, 'rank')
    values = double(info.radar.rank);
elseif isfield(info, 'packet_metadata') && isfield(info.packet_metadata, 'values') && ...
        isfield(info.packet_metadata.values, 'Rank') && isfield(info.packet_metadata.values, 'Signal_Type')
    p = info.packet_metadata.values;
    selected = p.Signal_Type(:) == 0;
    if isfield(p, 'selected_row_zero_based')
        selected = selected & p.selected_row_zero_based(:) >= 0 & p.selected_row_zero_based(:) < Na;
    end
    values = double(p.Rank(selected));
else
    error('echo_SWST:MissingRank', 'Rank is required to calculate the absolute RAW range axis.');
end
values = unique(values(:));
assert(numel(values) == 1 && isfinite(values) && values >= 0 && fix(values) == values, ...
    'echo_SWST:VariableRank', 'Selected ECHO rows must share one finite nonnegative integer Rank.');
rank = values;
end

function removeTemporary(path)
if isfile(path), delete(path); end
end
