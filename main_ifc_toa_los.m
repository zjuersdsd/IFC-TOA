%% IFC-TOA: calibration and TOA estimation for the shared 32.122 m LoS recording
% Run this complete script. All paths are relative to this release folder.
% See README.md for parameter definitions, time origins and output units.
clear;
clc;

%% 1. Recording, acquisition and synchronization parameters
releaseDir = fileparts(mfilename('fullpath'));
addpath(fullfile(releaseDir, 'src'), '-begin');
audioFile = fullfile(releaseDir, 'data', 'los_32.122m.wav');
outputDir = fullfile(releaseDir, 'results');
fs = 48000;
temperatureC = 25;
soundSpeed = 331.3 + 0.606 * temperatureC;
distanceCalib = 0.679;                  % Known calibration distance, m
distanceTrue = 32.122;                  % Evaluation only; never used for TOA selection
recordingStart = 0.165;                 % Initial crop in the original WAV, s
calibTimeRange = [0, 20];               % Calibration interval after the initial crop, s
testTimeRange = [50, 80];               % Nominal measurement interval, s
testStart = 50.14 - recordingStart + 0.2; % Exact measurement start after cropping, s
framePeriod = 0.2;                     % Repetition period, not chirp duration, s
calibrationPeakThreshold = 0.3;        % Relative peak height for synchronization only
maxFrames = Inf;                       % Inf processes all 150 measurement frames

%% 2. Chirp and common preprocessing parameters
RefParameter.fs = fs;
RefParameter.dura = 0.03;              % T, s
RefParameter.B = 3000;                 % B, Hz
RefParameter.f0 = 19000;               % f0, Hz
RefParameter.f1 = RefParameter.f0 + RefParameter.B;
RefParameter.len = round(RefParameter.dura * fs);
RefParameter.k = RefParameter.B / RefParameter.dura;
bandpassOrder = 8;
bandpassHalfPowerHz = [18000, 23000];

%% 3. IFC-TOA parameters, unchanged from the fast experimental implementation
aif_opts.globalMaxLagSamples = round(0.2 * fs);
aif_opts.cropPreSamples = round(0.001 * fs); % tau_a = 1 ms
aif_opts.cropPostSamples = round(0.002 * fs);
aif_opts.localMaxLagSamples = round(0.033 * fs);
aif_opts.fOffset = 0;
aif_opts.lowpassCutoff = 5000;
aif_opts.voicesPerOctave = 48;
aif_opts.rho = 1;
aif_opts.w = round(0.001 * fs);             % N_w = 48 samples
aif_opts.beta = 0.1;                       % Relative GCC candidate threshold
aif_opts.dfTol = 20;                       % f_tol, Hz
aif_opts.dfTolMin = 20;
aif_opts.minRidgeBins = 3;
aif_opts.minValidWindowRatio = 0.5;
aif_opts.fbMin = 0;
aif_opts.fbPad = 500;
aif_opts.energyMetricThr = 0.15;            % eta_IFC
aif_opts.isNormal = false;
aif_opts.useStabilityWeight = true;

%% 4. Generate the reference and preprocess the complete recording
validateattributes(maxFrames, {'numeric'}, {'scalar', 'positive'});
assert(isinf(maxFrames) || (isfinite(maxFrames) && maxFrames == floor(maxFrames)), ...
    'maxFrames must be a positive integer or Inf.');
frameLen = round(framePeriod * fs);
nCalib = floor(diff(calibTimeRange) / framePeriod);
nFrames = min(floor(diff(testTimeRange) / framePeriod), maxFrames);
t = (0:RefParameter.len-1) / fs;
ref = cos(2*pi .* (RefParameter.f0 + RefParameter.k/2 .* t) .* t).';
[audio, fsRead] = audioread(audioFile);
assert(fsRead == fs && size(audio, 2) == 1, 'Expected a mono 48 kHz WAV.');
% Preserve the experimental indexing convention: no +1 for this initial crop.
audio = audio(round(recordingStart * fs):end);
bpFilt = designfilt('bandpassiir', 'FilterOrder', bandpassOrder, ...
    'HalfPowerFrequency1', bandpassHalfPowerHz(1), ...
    'HalfPowerFrequency2', bandpassHalfPowerHz(2), 'SampleRate', fs);
audio = filtfilt(bpFilt, audio);
lastTestSample = round(testStart * fs) + nFrames * frameLen;
assert(lastTestSample <= numel(audio), 'Insufficient recording length.');

%% 5. Estimate linear clock drift and the calibration timing reference
calibAudio = audio(round(calibTimeRange(1)*fs)+1:round(calibTimeRange(2)*fs));
calibrationLags = zeros(nCalib, 1);
for i = 1:nCalib
    sig = calibAudio((i-1)*frameLen+1:i*frameLen);
    calibrationLags(i) = calibration_arrival(sig, ref, calibrationPeakThreshold);
end
fitCoefficients = polyfit(1:nCalib, calibrationLags, 1);
driftSamplesPerFrame = fitCoefficients(1);
correctedCalibrationLags = calibrationLags - (0:nCalib-1).' * driftSamplesPerFrame;
syncReferenceSamples = mean(correctedCalibrationLags);
testOffsetSeconds = testStart - testTimeRange(1);

%% 6. Estimate TOA using only the fast IFC-TOA implementation
testAudio = audio(round(testStart*fs)+1:lastTestSample);
rawTOASamples = zeros(nFrames, 1);
for i = 1:nFrames
    sig = testAudio((i-1)*frameLen+1:i*frameLen);
    rawTOASamples(i) = aif_gcc_toa_estimator_real_fast(sig, ref, RefParameter, aif_opts);
    if mod(i, 25) == 0 || i == nFrames
        fprintf('IFC-TOA: %d/%d frames processed.\n', i, nFrames);
    end
end
frameIndex = (1:nFrames).';
rawTOASeconds = rawTOASamples / fs;
% Keep the original experiment's offset and nominal frame-count convention.
alignedSamples = abs(testOffsetSeconds * fs) + rawTOASamples - ...
    ((frameIndex-1) + testTimeRange(1)/framePeriod) * driftSamplesPerFrame;
relativeTOASeconds = (alignedSamples - syncReferenceSamples) / fs;
propagationTOASeconds = relativeTOASeconds + distanceCalib / soundSpeed;
estimatedDistance = relativeTOASeconds * soundSpeed + distanceCalib;
absoluteError = abs(estimatedDistance - distanceTrue);

%% 7. Save TOA and ranging results
resultTable = table(frameIndex, rawTOASamples, rawTOASeconds, relativeTOASeconds, ...
    propagationTOASeconds, estimatedDistance, absoluteError, ...
    'VariableNames', {'Frame', 'RawTOA_samples', 'RawTOA_s', 'RelativeTOA_s', ...
    'PropagationTOA_s', 'EstimatedDistance_m', 'AbsoluteError_m'});
summaryTable = table(distanceTrue, nFrames, mean(absoluteError), median(absoluteError), ...
    prctile(absoluteError, 95), 'VariableNames', ...
    {'GroundTruth_m', 'Frames', 'MAE_m', 'MedianAbsoluteError_m', 'P95AbsoluteError_m'});
settings = struct('fs', fs, 'temperatureC', temperatureC, 'soundSpeed', soundSpeed, ...
    'distanceCalib', distanceCalib, 'distanceTrue', distanceTrue, ...
    'recordingStart', recordingStart, 'calibTimeRange', calibTimeRange, ...
    'testTimeRange', testTimeRange, 'testStart', testStart, 'framePeriod', framePeriod, ...
    'calibrationPeakThreshold', calibrationPeakThreshold, 'maxFrames', maxFrames, ...
    'bandpassOrder', bandpassOrder, 'bandpassHalfPowerHz', bandpassHalfPowerHz);
if ~isfolder(outputDir)
    mkdir(outputDir);
end
writetable(resultTable, fullfile(outputDir, 'toa_per_frame.csv'));
writetable(summaryTable, fullfile(outputDir, 'summary.csv'));
save(fullfile(outputDir, 'ifc_toa_results.mat'), 'resultTable', 'summaryTable', ...
    'settings', 'RefParameter', 'aif_opts', 'calibrationLags', ...
    'correctedCalibrationLags', 'fitCoefficients', 'driftSamplesPerFrame', ...
    'syncReferenceSamples', 'testOffsetSeconds');
disp(summaryTable);
fprintf('Results saved to: %s\n', outputDir);

function lagSample = calibration_arrival(sig, ref, relativeThreshold)
% Synchronization helper: retain the original calibration selection exactly.
% Calibration searches the full correlation lag axis, unlike IFC candidate search.
    [cc, lags] = xcorr(sig, ref);
    magnitude = abs(cc);
    [~, locations] = findpeaks(magnitude, 'MinPeakHeight', relativeThreshold * max(magnitude));
    assert(~isempty(locations), 'No calibration peak found. Check recording alignment.');
    lagSample = lags(locations(1));
end
