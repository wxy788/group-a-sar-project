function [echo, info] = read_s1_raw(inputFile, numPulses, varargin)
%READ_S1_RAW Decode Sentinel-1 L0 RAW packets and load complex echoes in MATLAB.
%   [echo, info] = read_s1_raw() reads the HH image PRI span of the LA scene.
%   [echo, info] = read_s1_raw('', Inf, 'Polarization', 'HV');
%   [data, info] = read_s1_raw('', Inf, 'LoadIntoMemory', false);
%   block = data.iq(1:1000, :);  % Read part of the complete on-disk matrix.
%   [echo, info] = read_s1_raw('D:\data\measurement.dat', 2200, ...
%       'StartPulse', 1000, 'CarrierFrequencyHz', 5.405e9, 'Plot', false);
%   numPulses counts PRI slots including calibration/noise/missing zero rows.
%   StartPulse is one-based relative to the FIRST ECHO packet's PRI (default 1).
%   Inf reads through the LAST ECHO packet's PRI; finite counts stop there too.
%   Leading/trailing non-imaging slots are outside the image selection span.
%   info.image.Na_all is the full first-to-last ECHO PRI span, including gaps,
%   independent of StartPulse and numPulses. It is not the ECHO packet count.
%   CarrierFrequencyHz defaults to 5.405e9; info.radar.fc is Hz and
%   info.radar.lambda = 299792458 / fc is metres, also saved in the output MAT.
%   All ECHO chunks in the selected PRI span are decoded. Other slots are zero.
%   echo: complex single, rows=PRI slots, columns=relative range samples.
%   Its width is the largest 2*Number_of_Quads among ECHO packets in the whole
%   file, independent of the selection. Shorter echoes are zero-padded at the
%   end without scaling or resampling; each retains its original SWST origin.
%   NQ/SWL may vary; changes in sampling rate, swath or PRI are rejected.
%   info.timeline.echo_num_samples records each selected PRI's original valid
%   complex sample count (2*NQ for ECHO, zero for non-ECHO/missing slots).
%   echo_SWST retains this original length; after shifting it is not the last
%   valid column index. Use echo_SWST separately to align window origins.
%   With LoadIntoMemory=false, echo is a matfile object exposing variable iq.
%   Full local image: 56435 x 24424 complex single, approximately 10.27 GiB.
%   info.orbit: GPS time and ECEF position/velocity state vectors from RAW headers.
%   info.attitude: GPS time, scalar-first ECI/J2000 quaternion and body angular rate.
%   Both retain the full file's unique native navigation records.
%   info.orbit.pulse_position_ecef_m / pulse_velocity_ecef_mps are Npri x 3,
%   evaluated at every selected PRI, including calibration/noise/missing slots.
%   Position/velocity use cubic Hermite interpolation. Short endpoint constant-
%   velocity extrapolation is limited to one native interval and marked by
%   info.orbit.pulse_extrapolated_mask. Missing navigation yields NaN arrays.
%   These per-PRI fields are returned in info only; the decoded MAT is not rewritten.
%   time_relative_s uses the first source packet as origin; time_gps_s is NOT UTC.
%   Uses the external Python decoder; no MATLAB Python bridge is required.
%   This decodes RAW I/Q only; it does not focus or radiometrically calibrate.

if nargin < 1 || isempty(inputFile), inputFile = ''; end
if nargin < 2 || isempty(numPulses), numPulses = Inf; end
validateattributes(numPulses, {'numeric'}, {'scalar','positive','real','nonnan'});
assert(isinf(numPulses) || fix(numPulses) == numPulses, 'numPulses must be a positive integer or Inf.');
toolDir = fileparts(mfilename('fullpath'));
parser = inputParser;
addParameter(parser, 'Polarization', 'HH');
addParameter(parser, 'StartPulse', 1);
addParameter(parser, 'CarrierFrequencyHz', 5.405e9);
addParameter(parser, 'Chunk', []);
addParameter(parser, 'Python', fullfile(toolDir, '.venv', 'Scripts', 'python.exe'));
addParameter(parser, 'OutputFile', '');
addParameter(parser, 'Plot', true);
addParameter(parser, 'LoadIntoMemory', true, @(x)islogical(x) && isscalar(x));
parse(parser, varargin{:});
opt = parser.Results;
validateattributes(opt.StartPulse, {'numeric'}, {'scalar','real','integer','positive','finite','<=',flintmax}, ...
    mfilename, 'StartPulse');
validateattributes(opt.CarrierFrequencyHz, {'numeric'}, {'scalar','real','positive','finite'}, ...
    mfilename, 'CarrierFrequencyHz');
if ~isempty(opt.Chunk)
    error('Chunk selection is no longer supported. All ECHO chunks are placed on the complete PRI timeline.');
end
pol = upper(char(opt.Polarization));
assert(ismember(pol, {'HH','HV','VV','VH'}), 'Use HH, HV, VV or VH.');

if strlength(string(inputFile)) == 0
    dataRoot = fullfile(fileparts(toolDir), '哨兵1号数据', '洛杉矶');
    pattern = ['s1a-s2-raw-s-' lower(pol) ...
        '-20150109t014106-20150109t014141-004088-004f1f.dat'];
    hits = dir(fullfile(dataRoot, '**', pattern));
    assert(numel(hits) == 1, 'Expected exactly one local measurement file. Pass inputFile explicitly.');
    inputFile = fullfile(hits(1).folder, hits(1).name);
end
inputFile = char(inputFile);
assert(isfile(inputFile), 'Measurement file not found: %s', inputFile);
assert(isfile(opt.Python), 'Python environment not found. See README.md or pass ''Python'', executable.');

outputFile = char(opt.OutputFile);
if isempty(outputFile)
    [~, stem] = fileparts(inputFile);
    countText = 'all';
    if isfinite(numPulses), countText = sprintf('%d', numPulses); end
    outputFile = fullfile(toolDir, 'outputs', sprintf('%s_image_start%d_n%s_fc%.17g.mat', ...
        stem, opt.StartPulse, countText, double(opt.CarrierFrequencyHz)));
end
helper = fullfile(toolDir, 'decode_s1_raw.py');
cmd = sprintf('%s -X utf8 %s --input %s --output %s --start-pulse %d --carrier-frequency-hz %.17g', ...
    shellQuote(opt.Python), shellQuote(helper), shellQuote(inputFile), ...
    shellQuote(outputFile), opt.StartPulse, double(opt.CarrierFrequencyHz));
if isfinite(numPulses), cmd = sprintf('%s --pulses %d', cmd, numPulses); end
fprintf('Decoding ECHO chunks from image PRI %d (one-based; first ECHO is 1)...\n', opt.StartPulse);
[status, commandOutput] = system(cmd);
fprintf('%s\n', commandOutput);
assert(status == 0, 'RAW decoding failed. Read the decoder message above.');

info = load(outputFile, 'radar', 'packet_metadata', 'chunk_summary', 'processing_info', ...
    'timeline', 'orbit', 'attitude', 'image');
assert(info.processing_info.complete, 'The output file does not contain a completed decode.');
info.mat_file = outputFile;
info = interpolate_s1_orbit(info);
data = matfile(outputFile);
dimensions = size(data, 'iq');
assert(numel(dimensions) == 2, 'Expected a two-dimensional iq matrix.');
fprintf('Complete matrix: %d PRI slots x %d range samples (%.2f GiB complex single).\n', ...
    dimensions(1), dimensions(2), prod(double(dimensions))*8/2^30);
if opt.LoadIntoMemory
    fprintf('Loading matrix into MATLAB memory; use LoadIntoMemory=false for disk access.\n');
    echo = data.iq;
    assert(isa(echo, 'single'), 'Expected a single-precision matrix.');
else
    echo = data;
end

if opt.Plot
    % Decimate BEFORE abs/log10 so a full-scene preview never duplicates 10+ GiB.
    displayRows = 1:max(1, ceil(dimensions(1)/1500)):dimensions(1);
    displayColumns = 1:max(1, ceil(dimensions(2)/2000)):dimensions(2);
    if opt.LoadIntoMemory
        displayEcho = echo(displayRows, displayColumns);
    else
        displayEcho = data.iq(displayRows, displayColumns);
    end
    figure('Name', 'Sentinel-1 decoded RAW echo');
    displayTime = info.timeline.time_relative_s(displayRows);
    imagesc(displayColumns, displayTime, 20 * log10(max(abs(displayEcho), single(1e-12))));
    % Robust display stretch: zero samples must not dominate the colour scale.
    sample = abs(displayEcho(1:max(1, ceil(numel(displayEcho)/200000)):end));
    sample = sort(double(sample(sample > 0 & isfinite(sample))));
    if ~isempty(sample)
        limits = 20 * log10(sample(max(1, round([0.02, 0.98]*numel(sample)))));
        if limits(2) > limits(1), clim(limits); end
    end
    axis xy;
    colormap gray;
    bar = colorbar;
    bar.Label.String = 'Amplitude (dB, arbitrary reference)';
    xlabel('Range sample (fast time)');
    ylabel('Slow time from first file PRI (s)');
    title('RAW I/Q on PRI timeline (decimated display; zero slots retained)');
end
end

function quoted = shellQuote(value)
% Windows cmd expands percent expressions even inside double quotes.
value = char(value);
assert(~any(ismember(value, [char(34), char(10), char(13), '%', '!'])), ...
    'Paths containing quotes, percent signs, exclamation marks or newlines are unsupported.');
quoted = [char(34) value char(34)];
end
