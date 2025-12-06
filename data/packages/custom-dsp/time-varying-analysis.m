% time-varying-analysis.m
% Functions for time-varying spectral analysis

function [freqs, amps] = analyze_static_spectrum(audio_file, window_size)
    % Analyze a single spectral snapshot from an audio file
    %
    % Parameters:
    %   audio_file - path to .wav file (string)
    %   window_size - FFT window size (optional, default 4096)
    %
    % Returns:
    %   freqs - array of peak frequencies in Hz
    %   amps - array of peak amplitudes in dBFS

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


function [times, freqs_matrix, amps_matrix, metadata] = analyze_spectral_evolution(audio_file, window_size, hop_size, min_amplitude)
    % Analyze time-varying spectrum across an audio file
    %
    % Parameters:
    %   audio_file - path to .wav file (string)
    %   window_size - FFT window size (optional, default 4096)
    %   hop_size - samples between frames (optional, default window_size/2)
    %   min_amplitude - minimum dBFS to include partial (optional, default 15)
    %
    % Returns:
    %   times - 1D array of time points in seconds (row vector)
    %   freqs_matrix - 2D array [num_frames × max_partials], zero-padded
    %   amps_matrix - 2D array [num_frames × max_partials], zero-padded
    %   metadata - struct with fields:
    %     .num_frames - actual number of frames analyzed
    %     .num_partials_per_frame - 1D array of partial counts per frame
    %     .sample_rate - sample rate in Hz
    %     .window_size - FFT window size used
    %     .hop_size - hop size used

    % Handle default parameters
    if nargin < 2
        window_size = 4096;
    end
    if nargin < 3
        hop_size = floor(window_size / 2);
    end
    if nargin < 4
        min_amplitude = 15;
    end

    % Load and prepare audio
    [y, fs] = audioread(audio_file);

    % Convert to mono
    if size(y, 2) > 1
        y = y(:, 1);
    end

    % High-pass filter at 80 Hz
    [b, a] = butter(4, 80/(fs/2), 'high');
    y = filter(b, a, y);

    % Calculate frame parameters
    num_samples = length(y);
    num_frames = floor((num_samples - window_size) / hop_size) + 1;

    % Pre-allocate arrays
    max_partials = 30;  % Estimate
    freqs_matrix = zeros(num_frames, max_partials);
    amps_matrix = zeros(num_frames, max_partials);
    times = zeros(1, num_frames);
    partials_per_frame = zeros(1, num_frames);

    % Analyze each frame
    for frame_idx = 1:num_frames
        % Extract window
        start_idx = (frame_idx - 1) * hop_size + 1;
        end_idx = start_idx + window_size - 1;

        if end_idx > num_samples
            break;
        end

        frame = y(start_idx:end_idx);
        times(frame_idx) = (start_idx - 1) / fs;

        % Analyze this frame
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
end


function [peak_freqs, peak_amps] = find_harmonics_frame(frame, fs, window_size, min_amplitude)
    % Analyze a single frame for harmonic peaks (internal helper function)
    %
    % Parameters:
    %   frame - audio samples for this frame (column vector)
    %   fs - sample rate in Hz
    %   window_size - FFT window size
    %   min_amplitude - minimum dBFS threshold (optional, default 15)
    %
    % Returns:
    %   peak_freqs - frequencies of detected peaks in Hz
    %   peak_amps - amplitudes of detected peaks in dBFS

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

printf('Time-varying analysis functions loaded!\n');
