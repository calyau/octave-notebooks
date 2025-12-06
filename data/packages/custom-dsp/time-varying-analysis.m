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
        min_amplitude = -60;
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
    if nargin < 4
        min_amplitude = -60;
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
    valid = peak_amps > min_amplitude;
    %printf('  [DEBUG] %d peaks above %.1f dBFS threshold\n', sum(valid), min_amplitude);

    peak_freqs = peak_freqs(valid);
    peak_amps = peak_amps(valid);

    if isempty(peak_freqs)
        return;
    end

    % Sort by amplitude and keep top partials
    [~, sort_idx] = sort(peak_amps, 'descend');
    max_peaks = min(20, length(peak_freqs));

    peak_freqs = peak_freqs(sort_idx(1:max_peaks));
    peak_amps = peak_amps(sort_idx(1:max_peaks));

    % Sort by frequency for output
    [peak_freqs, sort_idx] = sort(peak_freqs);
    peak_amps = peak_amps(sort_idx);

    %printf('  [DEBUG] Returning %d partials\n', length(peak_freqs));
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

    hold off;
    xlabel('Time (s)');
    ylabel('Ratio to Fundamental');
    if ~isempty(title_str)
        title([title_str ' - Harmonic Ratio Evolution (Inharmonicity)']);
    else
        title('Harmonic Ratio Evolution (Inharmonicity)');
    end
    if ~isempty(legend_labels)
        legend(legend_labels, 'location', 'best');
    end
    if ~isempty(time_range)
        xlim(time_range);
    end

    % Set y-range (either custom or auto with reference lines)
    if ~isempty(ratio_range)
        ylim(ratio_range);
    else
        ylim([1.5 max_partial + 0.5]);
    end

    grid on;

    % Add reference lines for perfect harmonics
    current_ylim = ylim;
    for i = 2:max_partial
        if i >= current_ylim(1) && i <= current_ylim(2)
            yline(i, '--', 'Color', [0.5 0.5 0.5], 'Alpha', 0.3);
        end
    end
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


function plot_spectral_richness(times, amps_matrix, threshold, title_str, time_range, richness_range)
    % Plot number of active partials over time
    %
    % Args:
    %   times: time vector (seconds)
    %   amps_matrix: 2D array [frames × partials] of amplitudes (dBFS)
    %   threshold: minimum dB to count as "active" (default -60)
    %   title_str: optional title prefix (default '')
    %   time_range: [min_time, max_time] to display (default: auto)
    %   richness_range: [min_count, max_count] of partials (default: auto)

    if nargin < 3
        threshold = -60;
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

    % Count active partials per frame
    active_partials = sum(amps_matrix > threshold, 2);

    figure;
    plot(times, active_partials, 'LineWidth', 2, 'Color', [0.8 0.2 0.2]);
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
        ylim([0 max(active_partials) + 2]);
    end
    grid on;

    % Add shading
    hold on;
    area(times, active_partials, 'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold off;
end


function spectral_analysis_report(times, freqs_matrix, amps_matrix, title_str, num_partials)
    % Generate a comprehensive spectral analysis report with all plots
    %
    % Args:
    %   times: time vector (seconds)
    %   freqs_matrix: 2D array [frames × partials] of frequencies
    %   amps_matrix: 2D array [frames × partials] of amplitudes (dBFS)
    %   title_str: optional title prefix (default '')
    %   num_partials: number of partials to plot in trajectories (default 10)
    %
    % Example:
    %   [times, freqs, amps, meta] = analyze_spectral_evolution('audio/oboe.wav');
    %   spectral_analysis_report(times, freqs, amps, 'Oboe A4', 10);

    if nargin < 4
        title_str = '';
    end
    if nargin < 5
        num_partials = 10;
    end
    min_amplitude = -60;

    printf('\nGenerating spectral analysis report...\n');

    % Plot 1: Partial trajectories
    printf('  - Partial trajectories\n');
    plot_partial_trajectories(times, freqs_matrix, amps_matrix, num_partials, title_str);

    % Plot 2: Harmonic ratios
    printf('  - Harmonic ratios\n');
    plot_harmonic_ratios(times, freqs_matrix, 6, title_str);

    % Plot 3: Spectral centroid
    printf('  - Spectral centroid\n');
    plot_spectral_centroid(times, freqs_matrix, amps_matrix, title_str);

    % Plot 4: Spectral richness
    printf('  - Spectral richness\n');
    plot_spectral_richness(times, amps_matrix, min_amplitude, title_str);

    % Print summary statistics
    printf('\nSpectral Analysis Summary:\n');
    printf('  Duration: %.2f seconds\n', times(end) - times(1));
    printf('  Number of frames: %d\n', length(times));
    printf('  Time resolution: %.3f seconds\n', mean(diff(times)));

    active_per_frame = sum(amps_matrix > min_amplitude, 2);
    printf('  Active partials: %.1f (mean), %d (max)\n', ...
           mean(active_per_frame), max(active_per_frame));

    % Calculate mean spectral centroid
    num_frames = length(times);
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
