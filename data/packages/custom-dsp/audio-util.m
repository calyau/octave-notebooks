% ============================================================================
% AUDIO AND SPECTRAL DATA CLEANING UTILITIES
% ============================================================================

% Detect and trim silence from audio signal
%
% This function analyzes the raw audio waveform to find where actual sound
% (including soft attacks and decays) begins and ends, trimming only the
% absolute silence at the beginning and end of the recording.
%
% Parameters:
%   y - audio signal (mono or stereo)
%   fs - sample rate in Hz
%   options - optional struct with:
%     .frame_size - analysis frame size in samples (default 2048)
%     .hop_size - hop between frames in samples (default 1024)
%     .min_duration - minimum duration of sound in seconds (default 0.1)
%     .silence_floor_db - dB value for completely silent frames (default -120)
%     .threshold_db - dB below peak to consider as silence (default -80)
%
% Returns:
%   y_trimmed - trimmed audio signal
%   fs - sample rate (unchanged)
%   trim_info - struct with trimming statistics

function [y_trimmed, fs, trim_info] = find_sound_boundaries(y, fs, options)
    % Handle default options
    if nargin < 3
        options = struct();
    end
    if ~isfield(options, 'frame_size')
        options.frame_size = 2048;
    end
    if ~isfield(options, 'hop_size')
        options.hop_size = 1024;
    end
    if ~isfield(options, 'min_duration')
        options.min_duration = 0.1;
    end
    if ~isfield(options, 'silence_floor_db')
        options.silence_floor_db = -120;
    end
    if ~isfield(options, 'threshold_db')
        options.threshold_db = -80;
    end

    % Convert to mono if stereo
    if size(y, 2) > 1
        y = y(:, 1);
    end

    original_length = length(y);

    % Calculate RMS energy in overlapping frames
    frame_size = options.frame_size;
    hop_size = options.hop_size;
    num_frames = floor((length(y) - frame_size) / hop_size) + 1;

    rms_db = zeros(num_frames, 1);
    frame_centers = zeros(num_frames, 1);

    for i = 1:num_frames
        start_idx = (i - 1) * hop_size + 1;
        end_idx = start_idx + frame_size - 1;

        if end_idx > length(y)
            break;
        end

        frame = y(start_idx:end_idx);
        rms_val = sqrt(mean(frame.^2));

        if rms_val > 0
            rms_db(i) = 20 * log10(rms_val);
        else
            rms_db(i) = options.silence_floor_db;
        end

        frame_centers(i) = floor((start_idx + end_idx) / 2);
    end

    % Trim to actual calculated frames
    actual_frames = sum(frame_centers > 0);
    rms_db = rms_db(1:actual_frames);
    frame_centers = frame_centers(1:actual_frames);

    % Find peak RMS and calculate threshold
    peak_db = max(rms_db);
    threshold = peak_db + options.threshold_db;

    % Find frames above threshold
    active_frames = find(rms_db > threshold);

    if isempty(active_frames)
        printf('Warning: No sound detected above threshold (%.1f dB)\n', threshold);
        y_trimmed = y;
        trim_info.start_sample = 1;
        trim_info.end_sample = length(y);
        trim_info.start_time = 0;
        trim_info.end_time = (length(y) - 1) / fs;
        trim_info.samples_removed_start = 0;
        trim_info.samples_removed_end = 0;
        trim_info.original_length = original_length;
        trim_info.trimmed_length = length(y);
        return;
    end

    % Find start and end of sound
    first_active_frame = active_frames(1);
    last_active_frame = active_frames(end);

    % Map back to sample indices
    start_sample = max(1, (first_active_frame - 1) * hop_size + 1);
    end_sample = min(length(y), last_active_frame * hop_size + frame_size);

    % Trim the signal
    y_trimmed = y(start_sample:end_sample);

    % Build info struct
    trim_info.start_sample = start_sample;
    trim_info.end_sample = end_sample;
    trim_info.start_time = (start_sample - 1) / fs;
    trim_info.end_time = (end_sample - 1) / fs;
    trim_info.samples_removed_start = start_sample - 1;
    trim_info.samples_removed_end = original_length - end_sample;
    trim_info.original_length = original_length;
    trim_info.trimmed_length = length(y_trimmed);
    trim_info.original_duration = (original_length - 1) / fs;
    trim_info.trimmed_duration = (length(y_trimmed) - 1) / fs;
    trim_info.threshold_used = threshold;
    trim_info.peak_db = peak_db;

    % Print summary
    printf('\n=== Audio Trimming Summary ===\n');
    printf('Original duration:  %.3f seconds (%d samples)\n', ...
           trim_info.original_duration, original_length);
    printf('Trimmed duration:   %.3f seconds (%d samples)\n', ...
           trim_info.trimmed_duration, trim_info.trimmed_length);
    printf('Sound starts at:    %.3f seconds (removed %.3f s from start)\n', ...
           trim_info.start_time, trim_info.start_time);
    printf('Sound ends at:      %.3f seconds (removed %.3f s from end)\n', ...
           trim_info.end_time, trim_info.original_duration - trim_info.end_time);
    printf('Peak RMS level:     %.1f dB\n', peak_db);
    printf('Threshold used:     %.1f dB (peak %.1f dB)\n', ...
           threshold, options.threshold_db);
    printf('==============================\n\n');
end


% Remove silent frames from spectral analysis results
%
% This function post-processes the output from analyze_spectral_evolution()
% to remove frames where no partials are active (all zeros or below threshold).
% It also resets the time axis to start at 0.
%
% Parameters:
%   times - time vector from analyze_spectral_evolution (row vector)
%   freqs_matrix - frequency matrix [frames × partials]
%   amps_matrix - amplitude matrix [frames × partials]
%   options - optional struct with:
%     .min_active_partials - minimum number of active partials (default 1)
%     .threshold_db - minimum dB to consider active (default -60)
%
% Returns:
%   times_trimmed - trimmed and zero-reset time vector
%   freqs_trimmed - trimmed frequency matrix
%   amps_trimmed - trimmed amplitude matrix
%   trim_info - struct with trimming statistics

function [times_trimmed, freqs_trimmed, amps_trimmed, trim_info] = ...
    trim_spectral_data(times, freqs_matrix, amps_matrix, options)

    % Handle defaults
    if nargin < 4
        options = struct();
    end
    if ~isfield(options, 'min_active_partials')
        options.min_active_partials = 1;
    end
    if ~isfield(options, 'threshold_db')
        options.threshold_db = -60;
    end

    % Find frames where at least min_active_partials are above threshold
    active_per_frame = sum((amps_matrix ~= 0) & (amps_matrix > options.threshold_db), 2);
    valid_frames = active_per_frame >= options.min_active_partials;

    % Find the first and last valid frame
    valid_indices = find(valid_frames);

    if isempty(valid_indices)
        times_trimmed = [];
        freqs_trimmed = [];
        amps_trimmed = [];
        trim_info.frames_removed_start = 0;
        trim_info.frames_removed_end = 0;
        trim_info.original_frames = length(times);
        trim_info.trimmed_frames = 0;
        printf('Warning: No valid frames found above threshold %.1f dB\n', options.threshold_db);
        return;
    end

    first_valid = valid_indices(1);
    last_valid = valid_indices(end);

    % Trim all arrays to only include valid frames
    times_trimmed = times(first_valid:last_valid);
    freqs_trimmed = freqs_matrix(first_valid:last_valid, :);
    amps_trimmed = amps_matrix(first_valid:last_valid, :);

    % Reset time to start at 0
    time_offset = times_trimmed(1);
    times_trimmed = times_trimmed - time_offset;

    % Build info struct
    trim_info.first_valid_frame = first_valid;
    trim_info.last_valid_frame = last_valid;
    trim_info.frames_removed_start = first_valid - 1;
    trim_info.frames_removed_end = length(times) - last_valid;
    trim_info.original_frames = length(times);
    trim_info.trimmed_frames = length(times_trimmed);
    trim_info.time_offset = time_offset;
    trim_info.original_start_time = times(1);
    trim_info.original_end_time = times(end);
    trim_info.trimmed_start_time = times_trimmed(1);
    trim_info.trimmed_end_time = times_trimmed(end);
    trim_info.original_duration = times(end) - times(1);
    trim_info.trimmed_duration = times_trimmed(end) - times_trimmed(1);

    % Print summary
    printf('\n=== Spectral Data Trimming Summary ===\n');
    printf('Original frames:    %d (%.3f to %.3f s, duration %.3f s)\n', ...
           trim_info.original_frames, trim_info.original_start_time, ...
           trim_info.original_end_time, trim_info.original_duration);
    printf('Trimmed frames:     %d (%.3f to %.3f s, duration %.3f s)\n', ...
           trim_info.trimmed_frames, trim_info.trimmed_start_time, ...
           trim_info.trimmed_end_time, trim_info.trimmed_duration);
    printf('Frames removed:     %d from start, %d from end\n', ...
           trim_info.frames_removed_start, trim_info.frames_removed_end);
    printf('Time offset:        %.3f s (subtracted to reset to t=0)\n', ...
           trim_info.time_offset);
    printf('======================================\n\n');
end

% Remove invalid data from individual partials
%
% This function post-processes amplitude/frequency matrices to handle the fact
% that different partials become active and fade out at different times.
% Instead of plotting zeros (which create visual spikes), it marks them as
% invalid so plotting functions can skip them or treat them as gaps.
%
% Parameters:
%   freqs_matrix - frequency matrix [frames × partials]
%   amps_matrix - amplitude matrix [frames × partials]
%   options - optional struct with:
%     .num_partials - max partials to print details for (default 30)
%     .replace_with - value to replace invalid data (default NaN)
%     .threshold_db - minimum dB to consider valid (default -60)
%
% Returns:
%   freqs_clean - cleaned frequency matrix (corresponding values replaced)
%   amps_clean - cleaned amplitude matrix (zeros replaced with NaN)
%   validity_mask - logical matrix [frames × partials] (true = valid data)

function [freqs_clean, amps_clean, validity_mask] = ...
    clean_partial_trajectories(freqs_matrix, amps_matrix, options)

    % Handle defaults
    if nargin < 3
        options = struct();
    end
    if ~isfield(options, 'num_partials')
        options.num_partials = 30;
    end
    if ~isfield(options, 'replace_with')
        options.replace_with = NaN;  % Use NaN so plotting functions skip these
    end
    if ~isfield(options, 'threshold_db')
        options.threshold_db = -60;
    end

    % Create validity mask
    % A point is valid if it's non-zero AND above threshold
    validity_mask = (amps_matrix ~= 0) & (amps_matrix > options.threshold_db);

    % Copy the matrices
    amps_clean = amps_matrix;
    freqs_clean = freqs_matrix;

    % For each partial, replace ALL invalid data with NaN (including interior gaps)
    orig_partial_count = size(amps_matrix, 2);

    printf('\n=== Cleaning Individual Partial Trajectories ===\n');

    for p = 1:orig_partial_count
        invalid_frames = ~validity_mask(:, p);

        % Replace all invalid frames with NaN
        amps_clean(invalid_frames, p) = options.replace_with;
        freqs_clean(invalid_frames, p) = options.replace_with;

        % Count valid and invalid
        num_valid = sum(validity_mask(:, p));
        num_invalid = sum(invalid_frames);

        % Find first and last valid for reporting
        valid_indices = find(validity_mask(:, p));

        if isempty(valid_indices)
            if p <= options.num_partials
                printf('Partial %2d: No valid data\n', p);
            end
            continue;
        end

        first_valid = valid_indices(1);
        last_valid = valid_indices(end);
        span = last_valid - first_valid + 1;
        gaps = span - num_valid;  % Interior gaps

        % Print summary for first N partials
        if p <= options.num_partials
            if gaps > 0
                printf('Partial %2d: %3d valid, %3d invalid (%d leading + %d trailing + %d gaps)\n', ...
                       p, num_valid, num_invalid, first_valid-1, ...
                       size(amps_matrix,1)-last_valid, gaps);
            else
                printf('Partial %2d: %3d valid, %3d invalid (%d leading + %d trailing)\n', ...
                       p, num_valid, num_invalid, first_valid-1, ...
                       size(amps_matrix,1)-last_valid);
            end
        end
    end

    printf('==============================================\n\n');
end

printf('Audio utility functions loaded ✓\n');
