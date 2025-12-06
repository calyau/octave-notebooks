# Implementation Instructions: Time-Varying Spectral Analysis

## Overview

Add time-varying spectral analysis capabilities to the existing Octave DSP package (`data/packages/custom-dsp/`). This will complement the existing point-in-time analysis with functions that track spectral evolution over time.

## Context

The existing `freq-analysis.m` file contains:
- `find_harmonics()` - analyzes a single spectral snapshot from the middle of an audio file
- `print_harmonics()` - displays harmonic analysis with MIDI information
- `plot_spectrum()` - plots frequency spectrum
- `plot_spectrogram()` - plots time-frequency representation
- Helper functions for MIDI conversion

## New Files to Create

### 1. `data/packages/custom-dsp/spectral-io.m`

Create a new file for saving/loading spectral data in .mat format (MATLAB v5/v6 compatible).

**Functions to implement:**

#### `save_static_spectrum(filename, freqs, amps, metadata)`
Save a single spectral snapshot to .mat file.

**Parameters:**
- `filename` - output .mat file path (string)
- `freqs` - 1D array of frequencies in Hz (row vector)
- `amps` - 1D array of amplitudes in dBFS (row vector)
- `metadata` - struct with fields:
  - `.fundamental` (required, scalar) - fundamental frequency in Hz
  - `.sample_rate` (required, scalar) - sample rate in Hz
  - `.instrument_id` (optional, scalar, default 0) - numeric instrument code
  - `.midi_note` (optional, scalar, default 0) - MIDI note number (0-127)
  - `.time_stamp` (optional, scalar, default 0.0) - time in seconds

**Implementation notes:**
- Reshape `freqs` and `amps` to row vectors using `reshape(x, 1, [])`
- Calculate `num_partials = length(freqs)`
- Use defaults for optional metadata fields with `isfield()` checks
- Save using: `save('-v6', filename, 'freqs', 'amps', 'num_partials', 'fundamental', 'sample_rate', 'time_stamp', 'instrument_id', 'midi_note')`
- Print confirmation: `printf('Saved static spectrum: %s (%d partials)\n', filename, num_partials)`

#### `save_spectral_evolution(filename, times, freqs_matrix, amps_matrix, metadata)`
Save time-varying spectral data to .mat file.

**Parameters:**
- `filename` - output .mat file path (string)
- `times` - 1D array of time points in seconds (row vector)
- `freqs_matrix` - 2D array [num_frames × max_partials], zero-padded
- `amps_matrix` - 2D array [num_frames × max_partials], zero-padded
- `metadata` - struct with fields:
  - `.num_partials_per_frame` (required, 1D array) - actual number of partials in each frame
  - `.fundamental` (required, scalar) - fundamental frequency in Hz
  - `.sample_rate` (required, scalar) - sample rate in Hz
  - `.window_size` (required, scalar) - FFT window size in samples
  - `.hop_size` (required, scalar) - hop size in samples between frames
  - `.instrument_id` (optional, scalar, default 0)
  - `.midi_note` (optional, scalar, default 0)

**Implementation notes:**
- Reshape `times` to row vector
- Reshape `num_partials_per_frame` to row vector
- Calculate `num_frames = length(times)`
- Use defaults for optional fields
- Rename matrices for saving: `freqs = freqs_matrix; amps = amps_matrix;`
- Save using: `save('-v6', filename, 'times', 'freqs', 'amps', 'num_frames', 'num_partials_per_frame', 'fundamental', 'sample_rate', 'window_size', 'hop_size', 'instrument_id', 'midi_note')`
- Print confirmation: `printf('Saved spectral evolution: %s (%d frames, %.2f-%.2f sec)\n', filename, num_frames, times(1), times(end))`

**Add at end of file:**
```octave
printf('Spectral I/O functions loaded!\n');
```

---

### 2. `data/packages/custom-dsp/time-varying-analysis.m`

Create a new file for time-varying spectral analysis functions.

**Functions to implement:**

#### `[freqs, amps] = analyze_static_spectrum(audio_file, window_size)`
Analyze a single spectral snapshot from an audio file.

**Parameters:**
- `audio_file` - path to .wav file (string)
- `window_size` - FFT window size (optional, default 4096)

**Returns:**
- `freqs` - array of peak frequencies in Hz
- `amps` - array of peak amplitudes in dBFS

**Implementation:**
```octave
function [freqs, amps] = analyze_static_spectrum(audio_file, window_size)
    if nargin < 2
        window_size = 4096;
    end
    
    [y, fs] = audioread(audio_file);
    
    % Convert to mono
    if size(y, 2) > 1
        y = y(:, 1);
    end
    
    % High-pass filter at 80 Hz
    [b, a] = butter(4, 80/(fs/2), 'high');
    y = filter(b, a, y);
    
    % Use existing find_harmonics function
    [freqs, amps] = find_harmonics(y, fs, window_size);
end
```

#### `[times, freqs_matrix, amps_matrix, metadata] = analyze_spectral_evolution(audio_file, window_size, hop_size, min_amplitude)`
Analyze time-varying spectrum across an audio file.

**Parameters:**
- `audio_file` - path to .wav file (string)
- `window_size` - FFT window size (optional, default 4096)
- `hop_size` - samples between frames (optional, default window_size/2)
- `min_amplitude` - minimum dBFS to include partial (optional, default 15)

**Returns:**
- `times` - 1D array of time points in seconds (row vector)
- `freqs_matrix` - 2D array [num_frames × max_partials], zero-padded
- `amps_matrix` - 2D array [num_frames × max_partials], zero-padded
- `metadata` - struct with fields:
  - `.num_frames` - actual number of frames analyzed
  - `.num_partials_per_frame` - 1D array of partial counts per frame
  - `.sample_rate` - sample rate in Hz
  - `.window_size` - FFT window size used
  - `.hop_size` - hop size used

**Implementation steps:**

1. **Load and prepare audio:**
   ```octave
   [y, fs] = audioread(audio_file);
   if size(y, 2) > 1
       y = y(:, 1);  % Mono
   end
   
   % High-pass filter
   [b, a] = butter(4, 80/(fs/2), 'high');
   y = filter(b, a, y);
   ```

2. **Calculate frame parameters:**
   ```octave
   num_samples = length(y);
   num_frames = floor((num_samples - window_size) / hop_size) + 1;
   ```

3. **Pre-allocate arrays:**
   ```octave
   max_partials = 30;  % Estimate
   freqs_matrix = zeros(num_frames, max_partials);
   amps_matrix = zeros(num_frames, max_partials);
   times = zeros(1, num_frames);
   partials_per_frame = zeros(1, num_frames);
   ```

4. **Analyze each frame:**
   ```octave
   for frame_idx = 1:num_frames
       % Extract window
       start_idx = (frame_idx - 1) * hop_size + 1;
       end_idx = start_idx + window_size - 1;
       
       if end_idx > num_samples
           break;
       end
       
       frame = y(start_idx:end_idx);
       times(frame_idx) = (start_idx - 1) / fs;
       
       % Analyze this frame (use helper function)
       [freqs, amps] = find_harmonics_frame(frame, fs, window_size, min_amplitude);
       
       num_found = length(freqs);
       partials_per_frame(frame_idx) = num_found;
       
       % Store (zero-padded)
       if num_found > 0
           freqs_matrix(frame_idx, 1:num_found) = freqs;
           amps_matrix(frame_idx, 1:num_found) = amps;
       end
       
       % Progress indicator every 10 frames
       if mod(frame_idx, 10) == 0
           printf('  Frame %d/%d (t=%.2fs)\n', frame_idx, num_frames, times(frame_idx));
       end
   end
   ```

5. **Trim and finalize:**
   ```octave
   % Trim if we broke early
   actual_frames = sum(times > 0);
   times = times(1:actual_frames);
   freqs_matrix = freqs_matrix(1:actual_frames, :);
   amps_matrix = amps_matrix(1:actual_frames, :);
   partials_per_frame = partials_per_frame(1:actual_frames);
   
   % Build metadata struct
   metadata.num_frames = actual_frames;
   metadata.num_partials_per_frame = partials_per_frame;
   metadata.sample_rate = fs;
   metadata.window_size = window_size;
   metadata.hop_size = hop_size;
   ```

#### `[peak_freqs, peak_amps] = find_harmonics_frame(frame, fs, window_size, min_amplitude)`
Analyze a single frame for harmonic peaks (internal helper function).

**Parameters:**
- `frame` - audio samples for this frame (column vector)
- `fs` - sample rate in Hz
- `window_size` - FFT window size
- `min_amplitude` - minimum dBFS threshold (optional, default 15)

**Returns:**
- `peak_freqs` - frequencies of detected peaks in Hz
- `peak_amps` - amplitudes of detected peaks in dBFS

**Implementation:**
```octave
function [peak_freqs, peak_amps] = find_harmonics_frame(frame, fs, window_size, min_amplitude)
    if nargin < 4
        min_amplitude = 15;
    end
    
    % Apply window
    windowed = frame .* hanning(length(frame));
    
    % Zero-pad
    padded = [windowed; zeros(window_size - length(frame), 1)];
    
    % Compute FFT and magnitude in dBFS
    spectrum = fft(padded);
    magnitude = abs(spectrum(1:floor(end/2)));
    magnitude = magnitude / (window_size / 2);  % Normalize for dBFS
    magnitude_db = 20 * log10(magnitude + 1e-10);
    
    freqs = (0:length(magnitude)-1) * fs / length(spectrum);
    
    % Find peaks above threshold
    peaks = [];
    for i = 2:length(magnitude_db)-1
        if magnitude_db(i) > magnitude_db(i-1) && ...
           magnitude_db(i) > magnitude_db(i+1) && ...
           magnitude_db(i) > min_amplitude
            peaks = [peaks; i];
        end
    end
    
    if isempty(peaks)
        peak_freqs = [];
        peak_amps = [];
        return;
    end
    
    % Sort by amplitude and keep top partials
    [sorted_amps, sort_idx] = sort(magnitude_db(peaks), 'descend');
    max_peaks = min(20, length(peaks));  % Keep top 20
    
    peak_indices = peaks(sort_idx(1:max_peaks));
    peak_freqs = freqs(peak_indices);
    peak_amps = magnitude_db(peak_indices);
    
    % Sort by frequency for output
    [peak_freqs, sort_idx] = sort(peak_freqs);
    peak_amps = peak_amps(sort_idx);
end
```

**Add at end of file:**
```octave
printf('Time-varying analysis functions loaded!\n');
```

---

## Usage Examples to Add

### 3. Create `data/packages/custom-dsp/examples/static-analysis-example.m`

```octave
% Example: Static spectrum analysis
pkg load signal;

% Load custom DSP functions
addpath('../');
freq_analysis;
spectral_io;
time_varying_analysis;

% Analyze oboe A4
[freqs, amps] = analyze_static_spectrum('../../audio/oboe.wav', 8192);

% Display results
printf('\nOboe A4 Static Analysis:\n');
print_harmonics(freqs, amps);

% Prepare metadata
meta.fundamental = 439.5;
meta.sample_rate = 48000;
meta.instrument_id = 1;  % 1 = oboe
meta.midi_note = 69;     % A4

% Save to .mat file
save_static_spectrum('../../output/oboe_A4_static.mat', freqs, amps, meta);
```

### 4. Create `data/packages/custom-dsp/examples/evolution-analysis-example.m`

```octave
% Example: Time-varying spectrum analysis
pkg load signal;

% Load custom DSP functions
addpath('../');
freq_analysis;
spectral_io;
time_varying_analysis;

% Analyze oboe A4 evolution
printf('Analyzing oboe evolution...\n');
[times, freqs, amps, meta] = analyze_spectral_evolution('../../audio/oboe.wav', 4096, 2048, 15);

% Add instrument metadata
meta.fundamental = 439.5;
meta.instrument_id = 1;  % oboe
meta.midi_note = 69;     % A4

% Save to .mat file
save_spectral_evolution('../../output/oboe_A4_evolution.mat', times, freqs, amps, meta);

printf('\nAnalysis complete!\n');
printf('Frames: %d\n', meta.num_frames);
printf('Duration: %.2f seconds\n', times(end) - times(1));
printf('Time resolution: %.3f seconds\n', meta.hop_size / meta.sample_rate);
```

---

## File Structure

After implementation, the structure should be:

```
data/
├── packages/
│   └── custom-dsp/
│       ├── freq-analysis.m          (existing)
│       ├── spectral-io.m            (NEW)
│       ├── time-varying-analysis.m  (NEW)
│       └── examples/
│           ├── static-analysis-example.m    (NEW)
│           └── evolution-analysis-example.m (NEW)
├── audio/
│   ├── oboe.wav
│   ├── clarinet.wav
│   └── together.wav
└── output/                          (create if doesn't exist)
```

---

## Testing Checklist

After implementation, test with:

1. **Static analysis:**
   ```octave
   cd data/packages/custom-dsp/examples
   static_analysis_example
   ```
   - Verify .mat file created in `data/output/`
   - Check file contains: `freqs`, `amps`, `num_partials`, `fundamental`, `sample_rate`, `time_stamp`, `instrument_id`, `midi_note`

2. **Evolution analysis:**
   ```octave
   evolution_analysis_example
   ```
   - Verify .mat file created
   - Check file contains: `times`, `freqs`, `amps`, `num_frames`, `num_partials_per_frame`, `fundamental`, `sample_rate`, `window_size`, `hop_size`, `instrument_id`, `midi_note`

3. **Verify .mat format:**
   ```octave
   whos -file ../../output/oboe_A4_static.mat
   whos -file ../../output/oboe_A4_evolution.mat
   ```
   - All variables should be type `double` (numeric arrays)
   - No strings or structs in the saved files

---

## Instrument ID Reference

For consistency across analyses:

```octave
% Instrument codes (for metadata.instrument_id):
% 0 = unspecified/unknown
% 1 = oboe
% 2 = clarinet
% 3 = flute
% 4 = piccolo
% 5 = violin
% 6 = viola
% 7 = cello
% 8 = synthesizer
```

---

## Notes

- All .mat files use MATLAB v5/v6 format (`-v6` flag) for maximum compatibility
- All saved data is numeric (no strings or structures in .mat files)
- Frequencies and amplitudes use consistent units: Hz and dBFS
- Zero-padding in evolution matrices ensures rectangular arrays
- Helper function `find_harmonics_frame()` is internal to `time-varying-analysis.m`
- Existing `find_harmonics()` function in `freq-analysis.m` is reused where appropriate

---

## Implementation Order

1. Create `spectral-io.m` with save functions
2. Create `time-varying-analysis.m` with analysis functions
3. Create example scripts in `examples/` directory
4. Create `output/` directory if it doesn't exist
5. Test static analysis
6. Test evolution analysis
7. Verify .mat file contents and format

---

## Expected Output Format

### Static .mat file variables:
- `freqs` - 1×N double (row vector)
- `amps` - 1×N double (row vector)
- `num_partials` - 1×1 double (scalar)
- `fundamental` - 1×1 double (scalar)
- `sample_rate` - 1×1 double (scalar)
- `time_stamp` - 1×1 double (scalar)
- `instrument_id` - 1×1 double (scalar)
- `midi_note` - 1×1 double (scalar)

### Evolution .mat file variables:
- `times` - 1×M double (row vector of M frames)
- `freqs` - M×K double (matrix: M frames × K max partials, zero-padded)
- `amps` - M×K double (matrix: M frames × K max partials, zero-padded)
- `num_frames` - 1×1 double (scalar)
- `num_partials_per_frame` - 1×M double (row vector)
- `fundamental` - 1×1 double (scalar)
- `sample_rate` - 1×1 double (scalar)
- `window_size` - 1×1 double (scalar)
- `hop_size` - 1×1 double (scalar)
- `instrument_id` - 1×1 double (scalar)
- `midi_note` - 1×1 double (scalar)
