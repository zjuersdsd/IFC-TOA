function [tauHatSample, out] = aif_gcc_toa_estimator_real_fast(y, s, params, opts)
% Fast IFC-TOA estimator for one preprocessed received frame.
% tauHatSample is a zero-based GCC lag in samples relative to the frame start.
% Use one output for early termination at the first accepted candidate.
% A second output evaluates all candidates and returns diagnostic arrays.
% The fixed analytic reference is cached; received signals are never cached.
% Parameters and processing conventions are documented in ../README.md.

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    if ~isfield(params, 'fs') || isempty(params.fs)
        error('params.fs must be specified.');
    end
    if ~isfield(params, 'k') || isempty(params.k)
        error('params.k must be specified.');
    end

    opts = set_default(opts, 'beta', 0.1);
    opts = set_default(opts, 'isNormal', false);

    Fs = params.fs;
    k  = params.k;

    y = y(:);
    s = s(:);

    Ny = length(y);
    Ns = length(s);

    opts = set_default(opts, 'globalMaxLagSamples', round(0.10 * Fs));
    opts = set_default(opts, 'cropPreSamples', round(0.005 * Fs));
    opts = set_default(opts, 'cropPostSamples', round(0.002 * Fs));
    opts = set_default(opts, 'localMaxLagSamples', opts.cropPreSamples + round(0.006 * Fs));


    opts = set_default(opts, 'fOffset', 0);
    opts = set_default(opts, 'lowpassCutoff', []);
    opts = set_default(opts, 'voicesPerOctave', 48);
    opts = set_default(opts, 'rho', 2);

    opts = set_default(opts, 'w', round(0.001 * Fs));
    opts = set_default(opts, 'minValidWindowRatio', 0.50);

    opts = set_default(opts, 'dfTol', []);
    opts = set_default(opts, 'dfTolMin', 40);
    opts = set_default(opts, 'minRidgeBins', 3);

    opts = set_default(opts, 'fbMin', 0);
    opts = set_default(opts, 'fbPad', 500);

    opts = set_default(opts, 'energyMetricThr', 0.25);
    opts = set_default(opts, 'useStabilityWeight', true);


    if isempty(opts.lowpassCutoff)
        if isfield(params, 'B') && ~isempty(params.B)
            opts.lowpassCutoff = min(0.45 * Fs, max(params.B, 3000));
        else
            opts.lowpassCutoff = min(0.45 * Fs, 5000);
        end
    end

    % Detect GCC candidates over nonnegative delays.
    [cc, lags] = xcorr(y, s);
    gccAbs = abs(cc);
    gccAbs = gccAbs / (max(gccAbs) + eps);

    posMask = lags >= 0 & lags <= opts.globalMaxLagSamples;
    lagPos = lags(posMask);
    gccPos = gccAbs(posMask);

    if isempty(lagPos)
        error('The nonnegative GCC search interval is empty.');
    end

    etaGcc = opts.beta*max(gccAbs);

    fpArgs = {'MinPeakHeight', etaGcc};

    [trustedPks, trustedLocs] = findpeaks(gccPos, fpArgs{:});
    trustedLags = lagPos(trustedLocs);

    if isempty(trustedLags)

        warning('No GCC candidate was found. Selecting the maximum GCC value.');
        [trustedPks, idxRough] = max(gccPos);
        trustedLags = lagPos(idxRough);
    end

    [roughLagSample, roughOrder] = min(trustedLags);
    roughPeakValue = trustedPks(roughOrder);

    % Advance cropping relative to the earliest GCC candidate.
    cropStartSample = max(0, roughLagSample - opts.cropPreSamples);

    cropLen = opts.localMaxLagSamples + opts.cropPostSamples;
    cropEndSample = min(Ny - 1, cropStartSample + cropLen - 1);

    cropIdx = (cropStartSample + 1) : (cropEndSample + 1);
    yCrop = y(cropIdx);
    segLen = length(yCrop);

    if segLen < max(16, round(0.5 * opts.w))
        error('The cropped signal is too short for SSWT analysis.');
    end

    if segLen < Ns
        warning('The cropped signal is shorter than the reference chirp.');
    end

    if isreal(yCrop)
        yCropA = hilbert(yCrop);
    else
        yCropA = yCrop;
    end

    % Cache only the analytic reference.
    persistent cachedReference cachedAnalyticReference
    if isempty(cachedReference) || ~isequal(cachedReference, s)
        cachedReference = s;
        if isreal(s)
            cachedAnalyticReference = hilbert(s);
        else
            cachedAnalyticReference = s;
        end
    end
    sA = cachedAnalyticReference;

    sPad = zeros(segLen, 1);
    Lcopy = min(Ns, segLen);
    sPad(1:Lcopy) = sA(1:Lcopy);

    tLocal = (0:segLen-1).' / Fs;

    % Dechirp; the theoretical beat frequency is k*n0/Fs + fOffset.
    u = conj(yCropA) .* sPad .* exp(1j * 2*pi*opts.fOffset*tLocal);

    fc = min(opts.lowpassCutoff, 0.45 * Fs);

    % Only the real part is needed for SSWT in the single-output path.
    u_lp = lowpass(real(u), fc, Fs);
    if nargout > 1
        u_lp = u_lp + 1j * lowpass(imag(u), fc, Fs);
    end

    % MATLAB R2022b default analytic Morlet wavelet, logarithmic scales.
    [sst_u, f_u] = wsst(real(u_lp), Fs, 'VoicesPerOctave', opts.voicesPerOctave);
    f_u = f_u(:);

    rho = opts.rho;
    mu = k / Fs;
    w = max(round(opts.w), 8);

    S = abs(sst_u).^rho;

    maxLocalLag = min(opts.localMaxLagSamples, segLen - 1);

    candidateLags = trustedLags(:);
    candidatePks = trustedPks(:);

    validCandMask = candidateLags >= cropStartSample & ...
                    candidateLags <= cropStartSample + maxLocalLag;

    candidateLags = candidateLags(validCandMask);
    candidatePks = candidatePks(validCandMask);

    if isempty(candidateLags)

        candidateLags = roughLagSample;
        candidatePks = roughPeakValue;
    end

    [candidateLags, sortIdx] = sort(candidateLags(:), 'ascend');
    candidatePks = candidatePks(sortIdx);

    nCandidate = candidateLags - cropStartSample;

    if isempty(opts.dfTol)
        df_sorted = diff(sort(f_u(:)));
        df_sorted = df_sorted(isfinite(df_sorted) & df_sorted > 0);
        if isempty(df_sorted)
            df0 = opts.dfTolMin;
        else
            df0 = median(df_sorted);
        end
        dfTol = max(3 * df0, opts.dfTolMin);
    else
        dfTol = opts.dfTol;
    end

    fbMin = opts.fbMin;
    fbMax = opts.fOffset + mu * maxLocalLag + opts.fbPad;

    beatMask = f_u >= fbMin & f_u <= fbMax;

    if ~any(beatMask)
        error('The beat analysis band contains no returned SSWT frequency bins.');
    end

    nCand = length(candidateLags);

    rawEnergyRatio = NaN(nCand, 1);
    energyMetric = NaN(nCand, 1);
    fbTheory = NaN(nCand, 1);
    ifMeanHz = NaN(nCand, 1);
    ifMadHz = NaN(nCand, 1);
    validWindowLength = zeros(nCand, 1);

    for jj = 1:nCand

        n0 = round(nCandidate(jj));

        if n0 < 0 || n0 > maxLocalLag
            continue;
        end

        timeIdx = (n0 + 1) : min(n0 + w, size(S, 2));
        validWindowLength(jj) = length(timeIdx);

        if length(timeIdx) < max(4, round(opts.minValidWindowRatio * w))
            continue;
        end

        rb = opts.fOffset + mu * n0;
        fbTheory(jj) = rb;

        ridgeMask = abs(f_u - rb) <= dfTol;

        if nnz(ridgeMask) < opts.minRidgeBins
            [~, idx0] = min(abs(f_u - rb));
            halfBins = floor(opts.minRidgeBins / 2);
            idx1 = max(1, idx0 - halfBins);
            idx2 = min(length(f_u), idx0 + halfBins);
            ridgeMask = false(size(f_u));
            ridgeMask(idx1:idx2) = true;
        end

        ridgeRows = find(ridgeMask);

        ridgeEnergy = sum(S(ridgeMask, timeIdx), 'all');
        totalEnergy = sum(S(beatMask, timeIdx), 'all') + eps;

        rawEnergyRatio(jj) = ridgeEnergy / totalEnergy;

        % Local IF: largest magnitude bin in each time column.
        [~, imax] = max(S(ridgeRows, timeIdx), [], 1);
        ifTrack = reshape(f_u(ridgeRows(imax)), 1, []);

        ifMeanHz(jj) = median(ifTrack, 'omitnan');
        ifMadHz(jj)  = mad(ifTrack, 1);

        freqBias = abs(ifMeanHz(jj) - rb);

        freqSpread = ifMadHz(jj);

        biasWeight = exp(-(freqBias / dfTol)^2);
        stableWeight = exp(-(freqSpread / dfTol)^2);

        if opts.useStabilityWeight
            energyMetric(jj) = rawEnergyRatio(jj) * biasWeight * stableWeight;
        else
            energyMetric(jj) = rawEnergyRatio(jj);
        end

        % Stop as soon as the earliest candidate passes IFC.
        if nargout < 2 && ~opts.isNormal && isfinite(energyMetric(jj)) && ...
                energyMetric(jj) >= opts.energyMetricThr
            tauHatSample = candidateLags(jj);
            return;
        end

    end
    if opts.isNormal
        energyMetric = energyMetric./(max(energyMetric)+eps);
    end

    etaPeak = opts.energyMetricThr;

    passIdx = find(isfinite(energyMetric) & energyMetric >= etaPeak, 1, 'first');

    if isempty(passIdx)

        warning('No candidate satisfies IFC. Selecting the earliest GCC candidate.');
        selectedIdx = 1;
    else
        selectedIdx = passIdx;
    end

    tauHatSample = candidateLags(selectedIdx);
    if nargout < 2
        return;
    end
    tauHatTime = tauHatSample / Fs;

    tauHatLocalSample = nCandidate(selectedIdx);
    tauHatLocalTime = tauHatLocalSample / Fs;

    tauMaxLocalSample = min(maxLocalLag, tauHatLocalSample + w - 1);
    tauMaxSample = cropStartSample + tauMaxLocalSample;
    tauMaxTime = tauMaxSample / Fs;

    if nargout > 1
        out = struct();

        out.tauHatSample = tauHatSample;
        out.tauHatTime = tauHatTime;

        out.tauHatLocalSample = tauHatLocalSample;
        out.tauHatLocalTime = tauHatLocalTime;

        out.tauMaxSample = tauMaxSample;
        out.tauMaxTime = tauMaxTime;
        out.tauMaxLocalSample = tauMaxLocalSample;
        out.tauMaxLocalTime = tauMaxLocalSample / Fs;

        out.roughLagSample = roughLagSample;
        out.roughLagTime = roughLagSample / Fs;
        out.roughPeakValue = roughPeakValue;

        out.cropStartSample = cropStartSample;
        out.cropStartTime = cropStartSample / Fs;
        out.cropEndSample = cropEndSample;
        out.cropEndTime = cropEndSample / Fs;

        out.gccAbs = gccAbs;
        out.lags = lags;
        out.gccPos = gccPos;
        out.lagPos = lagPos;

        out.etaGcc = etaGcc;
        out.trustedLags = trustedLags;
        out.trustedPks = trustedPks;

        out.candidateLags = candidateLags;
        out.candidateTimes = candidateLags / Fs;
        out.candidateLocalLags = nCandidate;
        out.candidatePks = candidatePks;

        out.rawEnergyRatio = rawEnergyRatio;
        out.energyMetric = energyMetric;
        out.etaPeak = etaPeak;

        out.fbTheory = fbTheory;
        out.ifMeanHz = ifMeanHz;
        out.ifMadHz = ifMadHz;
        out.validWindowLength = validWindowLength;

        out.selectedIdx = selectedIdx;

        out.yCrop = yCrop;
        out.u = u;
        out.u_lp = u_lp;
        out.sst_u = sst_u;
        out.f_u = f_u;

        out.dfTol = dfTol;
        out.fbMin = fbMin;
        out.fbMax = fbMax;

        out.params = opts;
        out.params.mu = mu;
        out.params.w = w;
    end
end

function opts = set_default(opts, fieldName, defaultValue)
    if ~isfield(opts, fieldName) || isempty(opts.(fieldName))
        opts.(fieldName) = defaultValue;
    end
end
