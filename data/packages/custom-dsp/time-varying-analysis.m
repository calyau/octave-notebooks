% time-varying-analysis.m
% Functions for time-varying spectral analysis

function [times, freqs_matrix, amps_matrix, metadata] = analyze_spectral_evolution(audio_file, options)
    % Analyze time-varying spectrum across an audio file
    %
    % Parameters:
    %   audio_file - path to .wav file (string)
    %   options - optional struct with:
    %     .high_pass_freq - high-pass filter cutoff in Hz (default 30)
    %     .hop_size - samples between frames (default window_size/2)
    %     .max_partials - maximum number of partials to detect per frame (default 30)
    %     .threshold_db - minimum dB to consider active (default: -60)
    %     .window_size - FFT window size (default 4096)
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
        options = struct();
    end

    if ~isfield(options, 'high_pass_freq')
        options.high_pass_freq = 20;
    end
    if ~isfield(options, 'max_partials')
        options.max_partials = 20;
    end
    if ~isfield(options, 'threshold_db')
        options.threshold_db = -60;
    end
    if ~isfield(options, 'window_size')
        options.window_size = 4096;
    end
    if ~isfield(options, 'hop_size')
        options.hop_size = floor(options.window_size / 2);
    end

    % Load and prepare audio
    [y, fs] = audioread(audio_file);

    % Convert to mono
    if size(y, 2) > 1
        y = y(:, 1);
    end

    % High-pass filter
    [b, a] = butter(4, options.high_pass_freq/(fs/2), 'high');
    y = filter(b, a, y);

    % Calculate frame parameters
    num_samples = length(y);
    num_frames = floor((num_samples - options.window_size) / options.hop_size) + 1;

    % Pre-allocate arrays
    freqs_matrix = zeros(num_frames, options.max_partials);
    amps_matrix = zeros(num_frames, options.max_partials);
    times = zeros(1, num_frames);
    partials_per_frame = zeros(1, num_frames);

    % Analyze each frame
    for frame_idx = 1:num_frames
        % Extract window
        start_idx = (frame_idx - 1) * options.hop_size + 1;
        end_idx = start_idx + options.window_size - 1;

        if end_idx > num_samples
            break;
        end

        frame = y(start_idx:end_idx);
        times(frame_idx) = (start_idx - 1) / fs;

        % Analyze this frame
        [freqs, amps] = find_harmonics_frame(frame, fs, options);

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
    metadata.window_size = options.window_size;
    metadata.hop_size = options.hop_size;

    % Detect fundamental using harmonic analysis
    metadata.fundamental = detect_fundamental(freqs_matrix, amps_matrix);
    if metadata.fundamental > 0
        metadata.midi_note = freq_to_midi_note(metadata.fundamental);
    else
        metadata.midi_note = 'N/A';
    end
end


function [peak_freqs, peak_amps] = find_harmonics_frame(frame, fs, options)
    % Handle default parameters
    if nargin < 3
        options = struct();
    end
    if ~isfield(options, 'max_partials')
        options.max_partials = 20;
    end
    if ~isfield(options, 'threshold_db')
        options.threshold_db = -60;
    end
    if ~isfield(options, 'window_size')
        options.window_size = 4096;
    end

    % Apply window
    windowed = frame .* hanning(length(frame));

    % Zero-pad
    padded = [windowed; zeros(options.window_size - length(frame), 1)];

    % Compute FFT and magnitude in dBFS
    spectrum = fft(padded);
    magnitude = abs(spectrum(1:floor(end/2)));
    magnitude = magnitude / (options.window_size / 2);  % Normalize for dBFS
    magnitude_db = 20 * log10(magnitude + 1e-10);

    freqs = (0:length(magnitude)-1)' * fs / length(spectrum);

    % DEBUG: Print some stats
    %printf('  [DEBUG] magnitude_db range: %.1f to %.1f dBFS\n', min(magnitude_db), max(magnitude_db));

    % Use findpeaks instead of manual peak detection
    % Shift spectrum to be all positive for findpeaks
    min_db = min(magnitude_db);
    spectrum_shifted = magnitude_db - min_db;

    % Calculate dynamic threshold
    peak_threshold = max(spectrum_shifted) - 42;
    if peak_threshold < 0
        peak_threshold = 0;
    end

    %printf('  [DEBUG] peak_threshold: %.1f (shifted), looking for peaks...\n', peak_threshold);

    % Find peaks
    [pks, locs] = findpeaks(spectrum_shifted, 'MinPeakHeight', peak_threshold, ...
                                              'MinPeakDistance', 10);

    %printf('  [DEBUG] Found %d peaks before amplitude filtering\n', length(locs));

    if isempty(locs)
        peak_freqs = [];
        peak_amps = [];
        return;
    end

    % Get actual dB values (not shifted)
    peak_freqs = freqs(locs);
    peak_amps = magnitude_db(locs);

    % Filter by minimum amplitude threshold
    valid = peak_amps > options.threshold_db;
    %printf('  [DEBUG] %d peaks above %.1f dBFS threshold\n', sum(valid), options.threshold_db);

    peak_freqs = peak_freqs(valid);
    peak_amps = peak_amps(valid);

    if isempty(peak_freqs)
        return;
    end

    % Sort by amplitude and keep top partials
    [~, sort_idx] = sort(peak_amps, 'descend');
    max_peaks = min(options.max_partials, length(peak_freqs));

    peak_freqs = peak_freqs(sort_idx(1:max_peaks));
    peak_amps = peak_amps(sort_idx(1:max_peaks));

    % Sort by frequency for output
    [peak_freqs, sort_idx] = sort(peak_freqs);
    peak_amps = peak_amps(sort_idx);

    %printf('  [DEBUG] Returning %d partials\n', length(peak_freqs));
end


function f0 = detect_fundamental(freqs_matrix, amps_matrix)
    % Detect fundamental frequency using harmonic analysis
    %
    % Algorithm:
    % 1. Skip first ~10% of frames (attack transient)
    % 2. Collect all detected peaks from stable frames
    % 3. For each candidate f0, score based on:
    %    - How many harmonics (n*f0) are present within tolerance
    %    - Weighted by amplitude of matching peaks
    % 4. Return the f0 with highest score

    num_frames = size(freqs_matrix, 1);
    if num_frames == 0
        f0 = 0;
        return;
    end

    % Use frames from stable portion (skip attack, use early-mid sustain)
    start_frame = max(1, floor(num_frames * 0.1));
    end_frame = min(num_frames, floor(num_frames * 0.5));
    if end_frame <= start_frame
        end_frame = num_frames;
    end

    % Collect all frequencies from stable frames
    stable_freqs = freqs_matrix(start_frame:end_frame, :);
    stable_amps = amps_matrix(start_frame:end_frame, :);
    all_freqs = stable_freqs(stable_freqs > 0);
    all_amps = stable_amps(stable_freqs > 0);

    if isempty(all_freqs)
        f0 = 0;
        return;
    end

    % Make sure they're column vectors
    all_freqs = all_freqs(:);
    all_amps = all_amps(:);

    % Get unique candidate fundamentals (all peaks are potential f0s)
    % Also consider sub-harmonics (peak/2, peak/3) as candidates
    candidates = unique([all_freqs; all_freqs/2; all_freqs/3]);
    candidates = candidates(candidates > 20 & candidates < 2000);  % Reasonable f0 range

    if isempty(candidates)
        f0 = 0;
        return;
    end

    best_score = 0;
    best_f0 = 0;

    for i = 1:length(candidates)
        f0_candidate = candidates(i);
        score = 0;

        % Check harmonics 1 through 10
        for harmonic = 1:10
            expected_freq = f0_candidate * harmonic;
            tolerance = expected_freq * 0.03;  % 3% tolerance

            % Find peaks near this expected harmonic
            matches = abs(all_freqs - expected_freq) < tolerance;
            if any(matches)
                % Add amplitude-weighted score (convert dB to linear)
                matching_amps = all_amps(matches);
                score = score + sum(10.^(matching_amps/20));
            end
        end

        if score > best_score
            best_score = score;
            best_f0 = f0_candidate;
        end
    end

    f0 = best_f0;
end


% ============================================================================
% VISUALIZATION FUNCTIONS FOR TIME-VARYING SPECTRAL ANALYSIS
% ============================================================================

function plot_partial_freq_trajectories(times, freqs_matrix, num_partials, title_str, time_range, freq_range)
    % Plot frequency trajectories for individual partials
    %
    % Args:
    %   times: time vector (seconds)
    %   freqs_matrix: 2D array [frames × partials] of frequencies
    %   num_partials: number of partials to plot (default 10)
    %   title_str: optional title prefix (default '')
    %   time_range: [min_time, max_time] to display (default: auto)
    %   freq_range: [min_freq, max_freq] for frequency plot (default: auto)

    if nargin < 3
        num_partials = 10;
    end
    if nargin < 4
        title_str = '';
    end
    if nargin < 5
        time_range = [];  % Auto range
    end
    if nargin < 6
        freq_range = [];  % Auto range
    end

    % Limit to available partials
    num_partials = min(num_partials, size(freqs_matrix, 2));

    % Create color map for partials
    colors = jet(num_partials);

    % Frequency trajectories
    figure;
    hold on;
    legend_labels = {};
    for partial_idx = 1:num_partials
        % Extract this partial's trajectory (skip zeros)
        trajectory = freqs_matrix(:, partial_idx);
        valid = trajectory > 0;  % Non-zero values

        if sum(valid) > 0
            plot(times(valid), trajectory(valid), 'o-', ...
                 'Color', colors(partial_idx,:), ...
                 'LineWidth', 1.5, ...
                 'MarkerSize', 4);
            legend_labels{end+1} = sprintf('Partial %d', partial_idx);
        end
    end
    hold off;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    if ~isempty(title_str)
        title([title_str ' - Partial Frequency Trajectories']);
    else
        title('Partial Frequency Trajectories');
    end
    if ~isempty(legend_labels)
        legend(legend_labels, 'location', 'eastoutside');
    end
    if ~isempty(time_range)
        xlim(time_range);
    end
    if ~isempty(freq_range)
        ylim(freq_range);
    end
    grid on;
end

function plot_partial_amp_trajectories(times, amps_matrix, num_partials, title_str, time_range, amp_range, margins, gaps)
    % Plot amplitude envelopes for individual partials, 4 per figure in 2x2 grid
    %
    % Args:
    %   times: time vector (seconds)
    %   amps_matrix: 2D array [frames × partials] of amplitudes
    %   num_partials: number of partials to plot (default 10)
    %   title_str: optional title prefix (default '')
    %   time_range: [min_time, max_time] to display (default: auto)
    %   amp_range: [min_amp, max_amp] for amplitude plot (default: auto)
    %   margins: [left, right, bottom, top] (default [0.08, 0.05, 0.08, 0.12])
    %   gaps: [horizontal, vertical] (default [0.08, 0.10])

    if nargin < 3
        num_partials = 10;
    end
    if nargin < 4
        title_str = '';
    end
    if nargin < 5
        time_range = [];
    end
    if nargin < 6
        amp_range = [];
    end
    if nargin < 7
        margins = [0.08, 0.05, 0.08, 0.12];
    end
    if nargin < 8
        gaps = [0.08, 0.10];
    end

    % Unpack margins and gaps
    left_margin = margins(1);
    right_margin = margins(2);
    bottom_margin = margins(3);
    top_margin = margins(4);

    horizontal_gap = gaps(1);
    vertical_gap = gaps(2);

    % Limit to available partials
    num_partials = min(num_partials, size(amps_matrix, 2));

    % Create color map for partials
    colors = jet(num_partials);

    % Plot 4 partials per figure (2x2 grid)
    partials_per_figure = 4;

    for start_idx = 1:partials_per_figure:num_partials
        figure;

        % Calculate which partials are in this figure
        end_idx = min(start_idx + partials_per_figure - 1, num_partials);

        % Create main title for this figure
        if ~isempty(title_str)
            if end_idx > start_idx
                suptitle_str = sprintf('%s - Partials %d-%d', title_str, start_idx, end_idx);
            else
                suptitle_str = sprintf('%s - Partial %d', title_str, start_idx);
            end
        else
            if end_idx > start_idx
                suptitle_str = sprintf('Partials %d-%d', start_idx, end_idx);
            else
                suptitle_str = sprintf('Partial %d', start_idx);
            end
        end

        annotation('textbox', [0 0.94 1 0.06], ...
                   'String', suptitle_str, ...
                   'EdgeColor', 'none', ...
                   'HorizontalAlignment', 'center', ...
                   'FontSize', 12, ...
                   'FontWeight', 'bold');

        % Calculate subplot dimensions
        total_width = 1 - left_margin - right_margin;
        total_height = 1 - top_margin - bottom_margin;
        subplot_width = (total_width - horizontal_gap) / 2;
        subplot_height = (total_height - vertical_gap) / 2;

        % Plot each partial in this figure (2x2 grid)
        for partial_idx = start_idx:end_idx
            % Calculate subplot position in 2x2 grid
            subplot_idx = partial_idx - start_idx + 1;

            % Calculate row and column (1-indexed)
            row = ceil(subplot_idx / 2);
            col = mod(subplot_idx - 1, 2) + 1;

            % Calculate position [left, bottom, width, height]
            left = left_margin + (col - 1) * (subplot_width + horizontal_gap);
            bottom = 1 - top_margin - (row * subplot_height) - ((row - 1) * vertical_gap);

            % Create subplot with custom position
            subplot('Position', [left, bottom, subplot_width, subplot_height]);

            plot_single_partial(times, amps_matrix, partial_idx, colors(partial_idx,:), ...
                               time_range, amp_range);
        end
    end
end

function plot_single_partial(times, amps_matrix, partial_idx, color, time_range, amp_range)
    % Helper function to plot a single partial

    trajectory = amps_matrix(:, partial_idx);
    valid = trajectory > -60;

    if sum(valid) > 0
        t_valid = times(valid);
        t_valid = t_valid(:);
        amp_valid = trajectory(valid);
        amp_valid = amp_valid(:);

        if ~isempty(amp_range)
            baseline = amp_range(1);
        else
            baseline = min(amp_valid) - 5;
        end

        baseline_vec = ones(size(amp_valid)) * baseline;

        hold on;
        % Fill area below curve
        fill([t_valid; flipud(t_valid)], ...
             [amp_valid; flipud(baseline_vec)], ...
             color, ...
             'FaceAlpha', 0.2, ...
             'EdgeColor', 'none');

        % Plot line
        plot(t_valid, amp_valid, '-', ...
             'Color', color, ...
             'LineWidth', 2);
        hold off;
    end

    if ~isempty(time_range)
        xlim(time_range);
    end
    if ~isempty(amp_range)
        ylim(amp_range);
    end

    xlabel('Time (s)');
    ylabel('Amplitude (dBFS)');
    title(sprintf('Partial %d', partial_idx));
    grid on;
end

function plot_harmonic_ratios(times, freqs_matrix, max_partial, title_str, time_range, ratio_range)
    % Plot harmonic ratios over time to show inharmonicity
    %
    % Args:
    %   times: time vector (seconds)
    %   freqs_matrix: 2D array [frames × partials] of frequencies
    %   max_partial: highest partial to plot (default 6)
    %   title_str: optional title prefix (default '')
    %   time_range: [min_time, max_time] to display (default: auto)
    %   ratio_range: [min_ratio, max_ratio] to display (default: auto)

    if nargin < 3
        max_partial = 6;
    end
    if nargin < 4
        title_str = '';
    end
    if nargin < 5
        time_range = [];
    end
    if nargin < 6
        ratio_range = [];
    end

    % Limit to available partials
    max_partial = min(max_partial, size(freqs_matrix, 2));

    figure;
    hold on;

    % Calculate harmonic ratios
    fundamental_trajectory = freqs_matrix(:, 1);

    colors = jet(max_partial - 1);

    legend_labels = {};
    for partial_idx = 2:max_partial
        ratios = freqs_matrix(:, partial_idx) ./ fundamental_trajectory;
        valid = ratios > 0 & ~isnan(ratios) & ~isinf(ratios);

        if sum(valid) > 0
            plot(times(valid), ratios(valid), 'o-', ...
                 'Color', colors(partial_idx-1,:), ...
                 'LineWidth', 1.5, ...
                 'MarkerSize', 4);
            legend_labels{end+1} = sprintf('Partial %d', partial_idx);
        end
    end

    % Set ranges
    if ~isempty(time_range)
        xlim(time_range);
    end

    if ~isempty(ratio_range)
        ylim(ratio_range);
    else
        ylim([1.5 max_partial + 0.5]);
    end

    current_ylim = ylim;
    current_xlim = xlim;

    % Add reference lines for perfect harmonics (using plot instead of yline)
    for i = 2:max_partial
        if i >= current_ylim(1) && i <= current_ylim(2)
            plot(current_xlim, [i i], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
        end
    end

    hold off;

    xlabel('Time (s)');
    ylabel('Ratio to Fundamental');

    % Multi-line title using newline
    if ~isempty(title_str)
        title(sprintf('%s\nHarmonic Ratio Evolution (Inharmonicity)', title_str));
    else
        title(sprintf('Harmonic Ratio Evolution\n(Inharmonicity)'));
    end

    if ~isempty(legend_labels)
        legend(legend_labels, 'location', 'eastoutside');
    end

    grid on;
end

function plot_spectral_centroid(times, freqs_matrix, amps_matrix, title_str, time_range, centroid_range)
    % Plot spectral centroid over time (shows "brightness")
    %
    % Args:
    %   times: time vector (seconds)
    %   freqs_matrix: 2D array [frames × partials] of frequencies
    %   amps_matrix: 2D array [frames × partials] of amplitudes (dBFS)
    %   title_str: optional title prefix (default '')
    %   time_range: [min_time, max_time] to display (default: auto)
    %   centroid_range: [min_centroid, max_centroid] in Hz (default: auto)

    if nargin < 4
        title_str = '';
    end
    if nargin < 5
        time_range = [];
    end
    if nargin < 6
        centroid_range = [];
    end

    % Calculate spectral centroid for each frame
    num_frames = length(times);
    centroids = zeros(num_frames, 1);

    for frame = 1:num_frames
        valid = freqs_matrix(frame, :) > 0;
        frame_freqs = freqs_matrix(frame, valid);
        frame_amps = amps_matrix(frame, valid);

        if ~isempty(frame_freqs)
            % Convert dB to linear amplitude
            linear_amps = 10 .^ (frame_amps / 20);

            % Weighted average frequency
            centroids(frame) = sum(frame_freqs .* linear_amps) / sum(linear_amps);
        end
    end

    figure;
    plot(times, centroids, 'LineWidth', 2, 'Color', [0.2 0.4 0.8]);
    xlabel('Time (s)');
    ylabel('Spectral Centroid (Hz)');
    if ~isempty(title_str)
        title([title_str ' - Timbral Brightness Over Time']);
    else
        title('Timbral Brightness Over Time');
    end
    if ~isempty(time_range)
        xlim(time_range);
    end
    if ~isempty(centroid_range)
        ylim(centroid_range);
    end
    grid on;

    % Add shading to show variation
    hold on;
    area(times, centroids, 'FaceColor', [0.2 0.4 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold off;
end


function plot_spectral_richness(times, amps_matrix, thresholds, title_str, time_range, richness_range)
    % Plot number of active partials over time with multiple threshold layers
    %
    % Args:
    %   times: time vector (seconds)
    %   amps_matrix: 2D array [frames × partials] of amplitudes (dBFS)
    %   thresholds: array of dB thresholds (default [-60, -50, -40, -30, -20])
    %   title_str: optional title prefix (default '')
    %   time_range: [min_time, max_time] to display (default: auto)
    %   richness_range: [min_count, max_count] of partials (default: auto)

    if nargin < 3 || isempty(thresholds)
        thresholds = [-60, -50, -40, -30, -20];
    end
    if nargin < 4
        title_str = '';
    end
    if nargin < 5
        time_range = [];
    end
    if nargin < 6
        richness_range = [];
    end

    % Ensure thresholds is a vector
    thresholds = thresholds(:)';  % Make it a row vector

    % Sort thresholds in ASCENDING order (lowest dB first = most partials)
    thresholds = sort(thresholds, 'ascend');

    num_thresholds = length(thresholds);

    % Create rainbow colors
    colors = jet(num_thresholds);

    % Calculate active partials for each threshold
    active_partials_all = zeros(length(times), num_thresholds);
    for i = 1:num_thresholds
        active_partials_all(:, i) = sum(amps_matrix > thresholds(i), 2);
    end

    % Find overall max for y-axis
    overall_max = max(active_partials_all(:));

    figure;
    hold on;

    % Plot areas from LOWEST threshold to HIGHEST (most partials to least)
    % This draws the largest area first, then smaller areas on top
    legend_labels = {};
    for i = 1:num_thresholds
        active_partials = active_partials_all(:, i);

        % Draw filled area
        area(times, active_partials, ...
             'FaceColor', colors(i,:), ...
             'FaceAlpha', 0.6, ...
             'EdgeColor', 'none');

        legend_labels{end+1} = sprintf('> %.0f dBFS', thresholds(i));
    end

    hold off;

    xlabel('Time (s)');
    ylabel('Number of Active Partials');
    if ~isempty(title_str)
        title([title_str ' - Spectral Richness Over Time']);
    else
        title('Spectral Richness Over Time');
    end
    if ~isempty(time_range)
        xlim(time_range);
    end
    if ~isempty(richness_range)
        ylim(richness_range);
    else
        ylim([0 overall_max + 2]);
    end

    % Add legend on the right
    legend(legend_labels, 'location', 'eastoutside');

    grid on;
end

function spectral_analysis_report(times, freqs_matrix, amps_matrix, title_str, num_partials)
    % Generate a comprehensive spectral analysis report with all plots
    %
    % Args:
    %   times: time vector (seconds)
    %   freqs_matrix: 2D array [frames × partials] of frequencies
    %   amps_matrix: 2D array [frames × partials] of amplitudes (dBFS)
    %   options - optional struct with:
    %     .title_prefix - title prefix (default '')
    %     .num_partials: number of partials to plot in trajectories (default all)
    %     .threshold_db - minimum dB to consider active (default: -60)
    %
    % Example:
    %   [times, freqs, amps, meta] = analyze_spectral_evolution('audio/oboe.wav');
    %   spectral_analysis_report(times, freqs, amps, 'Oboe A4', 10);
    frame_count = length(times);
    if nargin < 4
        options = struct();
    end
    if ~isfield(options, 'title_prefix')
        options.title_prefix = '';
    end
    if ~isfield(options, 'num_partials')
        options.num_partials = frame_count;
    end
    if options.num_partials == 'all'
        options.num_partials = frame_count;
    end
    if ~isfield(options, 'threshold_db')
        options.threshold_db = -60;
    end

    printf('\nGenerating spectral analysis report...\n');

    % Plot 1: Partial trajectories
    printf('  - Partial trajectories\n');
    plot_partial_trajectories(times, freqs_matrix, amps_matrix, options.num_partials, options.title_prefix);

    % Plot 2: Harmonic ratios
    printf('  - Harmonic ratios\n');
    plot_harmonic_ratios(times, freqs_matrix, 6, options.title_prefix);

    % Plot 3: Spectral centroid
    printf('  - Spectral centroid\n');
    plot_spectral_centroid(times, freqs_matrix, amps_matrix, options.title_prefix);

    % Plot 4: Spectral richness
    printf('  - Spectral richness\n');
    plot_spectral_richness(times, amps_matrix, options.threshold_db, options.title_prefix);

    % Print summary statistics
    printf('\nSpectral Analysis Summary:\n');
    printf('  Duration: %.2f seconds\n', times(end) - times(1));
    printf('  Number of frames: %d\n', frame_count);
    printf('  Time resolution: %.3f seconds\n', mean(diff(times)));

    active_per_frame = sum(amps_matrix > options.threshold_db, 2);
    printf('  Active partials: %.1f (mean), %d (max)\n', ...
           mean(active_per_frame), max(active_per_frame));

    % Calculate mean spectral centroid
    num_frames = frame_count;
    centroids = zeros(num_frames, 1);
    for frame = 1:num_frames
        valid = freqs_matrix(frame, :) > 0;
        frame_freqs = freqs_matrix(frame, valid);
        frame_amps = amps_matrix(frame, valid);
        if ~isempty(frame_freqs)
            linear_amps = 10 .^ (frame_amps / 20);
            centroids(frame) = sum(frame_freqs .* linear_amps) / sum(linear_amps);
        end
    end
    printf('  Spectral centroid: %.1f Hz (mean)\n', mean(centroids(centroids > 0)));

    printf('\nReport complete!\n');
end

printf('Time-varying analysis functions loaded!\n');
