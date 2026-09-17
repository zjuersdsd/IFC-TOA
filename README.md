# IFC-TOA: Public Code for a Single-Distance LoS Experiment

This folder provides a LoS recording and a complete processing workflow for the paper **Robust TOA Estimation for Acoustic Signals Based on Instantaneous Frequency Consistency in Multipath Environments**. The measurement distance is **32.122 m**, and the calibration distance is **0.679 m**. The main script uses the fast IFC-TOA implementation to process 150 measurement frames by default and outputs TOA estimates, propagation times after synchronization compensation, and ranging errors.

The parameter descriptions follow Sections III-B, V-A, and V-B of the revised manuscript. Recording alignment and auxiliary processing settings also follow the experimental code. This release contains the materials needed to reproduce the results for the shared recording. It does not include all experimental data, implementations of comparison algorithms, or runtime and complexity measurements.

## 1. Files and Usage

```text
IFC_TOA_public/
  main_ifc_toa_los.m                     Main script; parameters in Sections 1-3
  src/aif_gcc_toa_estimator_real_fast.m   Fast IFC-TOA and internal helper functions
  data/los_32.122m.wav                   Original mono recording, including calibration
  README.md                             Parameters and usage instructions
  results/                              Outputs; rerunning overwrites matching files
```

Use MATLAB R2022b with Signal Processing Toolbox, Wavelet Toolbox, and Statistics and Machine Learning Toolbox. The code uses `wsst`, `lowpass`, `designfilt`, `filtfilt`, `findpeaks`, `hilbert`, `mad`, and `prctile`. Other MATLAB versions have not been verified.

Place this folder at any location, open `main_ifc_toa_los.m` in MATLAB, and click **Run**. Paths are resolved relative to the main script, so the original experimental directory and other author code are not required. Use the default parameters for the first run. A complete run processes 100 calibration frames and 150 measurement frames. For a brief check of the workflow, set `maxFrames` to `3`; the resulting statistics will then describe only the first three measurement frames.

The shared WAV is an unchanged copy of the complete original recording, without resampling or truncation. The script applies the original bandpass preprocessing before extracting the calibration and measurement segments, preserving the filtering boundaries used in the experimental script.

## 2. Recording and Synchronization Parameters

| Main script variable | Default value | Description |
|---|---:|---|
| `audioFile` | `data/los_32.122m.wav` | Path to the shared recording. |
| `outputDir` | `results` | Directory for CSV and MAT outputs. |
| `fs` | 48000 Hz | Sampling rate $F_s$; it must match the WAV sampling rate. |
| `temperatureC` | 25 °C | Temperature used to calculate the sound speed. |
| `soundSpeed` | 346.45 m/s | Calculated as `331.3 + 0.606*temperatureC`. |
| `distanceCalib` | 0.679 m | Known distance during the initial calibration, used to recover absolute propagation time and distance. |
| `distanceTrue` | 32.122 m | Ground-truth straight-line distance measured using a laser rangefinder.|
| `recordingStart` | 0.165 s | Initial cropping position in the original WAV. |
| `calibTimeRange` | [0, 20] s | Calibration interval after the initial crop, containing 100 frames. |
| `testTimeRange` | [50, 80] s | Nominal measurement interval, used to determine the frame count and nominal frame indices for drift compensation. |
| `testStart` | 50.175 s | Actual measurement start after the initial crop, calculated as `50.14 - recordingStart + 0.2`. |
| `framePeriod` | 0.2 s | Chirp repetition period and frame duration. |
| `calibrationPeakThreshold` | 0.3 | Relative cross-correlation peak height threshold used only to establish the synchronization reference. This is separate from the IFC candidate threshold. |
| `maxFrames` | `Inf` | Processes all 150 measurement frames by default. A positive integer limits the number of processed frames. |

`recordingStart`, `testStart`, and the time intervals are alignment metadata for this recording, rather than general algorithm parameters transferable to other recordings. When replacing the recording, update the known calibration distance and segment start times accordingly.

Synchronization follows the initial timing reference and linear drift correction described in Section V-A. The earliest cross-correlation local peak exceeding the relative threshold is detected in each calibration frame. A linear fit of the peak delay against the frame index estimates the drift per frame. The mean delay after drift removal provides the synchronization reference. The calibration helper preserves the original search over the full correlation lag axis, including negative and positive delays; IFC candidate detection during measurement searches only nonnegative delays.

During measurement, the script applies the original correction based on the actual start offset and nominal frame index, then subtracts the calibration reference. This gives the propagation time difference relative to the calibration distance. Adding `distanceCalib/soundSpeed` yields the propagation time estimate. The raw GCC delay within a frame must not be interpreted directly as acoustic propagation time before synchronization compensation.

## 3. Chirp and Preprocessing Parameters

| Main script variable | Default value | Description |
|---|---:|---|
| `RefParameter.fs` | 48000 Hz | Sampling rate of the reference signal. |
| `RefParameter.dura` | 0.03 s | Chirp duration $T$. |
| `RefParameter.B` | 3000 Hz | Chirp bandwidth $B$. |
| `RefParameter.f0` | 19000 Hz | Initial chirp frequency $f_0$. |
| `RefParameter.f1` | 22000 Hz | Final chirp frequency, calculated as $f_0+B$. |
| `RefParameter.len` | 1440 samples | Number of chirp samples, calculated from $T*F_s$. |
| `RefParameter.k` | 100000 Hz/s | Chirp rate $k=B/T$, which maps a candidate delay to its theoretical beat frequency. |


## 4. IFC-TOA Parameters

The following fields belong to `aif_opts` in the main script. The principal manuscript parameters are `beta`, `cropPreSamples`, `w`, `dfTol`, `energyMetricThr`, and the SSWT and lowpass filter settings. The remaining fields specify implementation limits and discrete processing rules.

| Field | Default value | Description and relation to the manuscript |
|---|---:|---|
| `globalMaxLagSamples` | 9600 samples | Upper limit of the nonnegative GCC delay search, corresponding to 0.2 s. |
| `beta` | 0.1 | Relative GCC peak height threshold. The actual threshold is `beta*max(gccAbs)` and defines the candidate peak set. |
| `cropPreSamples` | 48 samples | Cropping advance, corresponding to the advance time $τ_a$=0.001 s. |
| `cropPostSamples` | 96 samples | Additional margin at the end of the cropped segment, corresponding to 0.002 s. |
| `lowpassCutoff` | 5000 Hz | Passband edge frequency supplied to `lowpass` after dechirping, rather than a specified 3 dB cutoff frequency. |
| `voicesPerOctave` | 48 | Number of SSWT scales per octave. MATLAB R2022b uses the default analytic Morlet mother wavelet with logarithmically spaced scales. |
| `rho` | 1 | Exponent applied to the SSWT magnitude: `abs(sst_u).^rho`. The default uses magnitude rather than squared magnitude. |
| `w` | 48 samples | Statistical window length N_w, corresponding to 1 ms. The window starts at the candidate delay. |
| `dfTol` | 20 Hz | Fixed frequency tolerance $f_tol$, used to define the neighborhood of the theoretical beat frequency and calculate the consistency weights. |
| `dfTolMin` | 20 Hz | Lower bound used only when automatic tolerance selection is enabled with `dfTol=[]`. It is inactive under the default fixed tolerance of 20 Hz. |
| `minRidgeBins` | 3 | If the frequency neighborhood contains too few bins, use the bin nearest the theoretical beat frequency and its neighbors. Fewer than three bins may be available at a frequency-axis boundary. |
| `minValidWindowRatio` | 0.5 | Minimum valid fraction of the statistical window. The default requires at least 24 valid time samples. |
| `energyMetricThr` | 0.15 | IFC acceptance threshold $η_{IFC}$. The earliest candidate meeting this threshold is selected. |
| `isNormal` | `false` | Disables normalization of candidate IFC scores by their maximum; the threshold of 0.15 is applied directly. |
| `useStabilityWeight` | `true` | Uses the complete IFC metric: energy concentration multiplied by the frequency bias and frequency fluctuation weights. |



## 5. Outputs

Each row of `results/toa_per_frame.csv` corresponds to one measurement:

| Column | Description |
|---|---|
| `Frame` | Measurement frame index, starting from 1. |
| `RawTOA_samples` | GCC delay selected by IFC, in samples. This is a zero-based lag relative to the frame start, not a MATLAB array index. |
| `RawTOA_s` | Raw delay within the frame, equal to `RawTOA_samples/fs`. |
| `RelativeTOA_s` | Propagation time difference after synchronization and drift compensation, relative to the calibration distance. |
| `PropagationTOA_s` | Propagation TOA estimate after adding the known calibration propagation time. |
| `EstimatedDistance_m` | Distance estimate calculated from the compensated time difference and calibration distance. |
| `AbsoluteError_m` | Absolute ranging error relative to the geometric ground-truth distance of 32.122 m. |

`results/summary.csv` reports the frame count, MAE, median absolute error, and 95th percentile absolute error. `results/ifc_toa_results.mat` also stores the result tables, parameters, calibration delays, fitted coefficients, and synchronization reference for inspection of the complete calculation.

A complete run with the default settings in MATLAB R2022b produced the following results. These outputs are included in the folder for comparison:

| Measurement frames | MAE | Median absolute error | 95th percentile absolute error |
|---:|---:|---:|---:|
| 150 | 0.0875587138 m | 0.0821893240 m | 0.1862097900 m |
