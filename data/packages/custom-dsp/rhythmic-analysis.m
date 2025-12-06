% ============================================================================
% RHYTHMIC ANALYSIS FUNCTIONS FOR AMPLITUDE ENVELOPES
% ============================================================================

function smoothed = smooth_envelope(envelope, window_size)
    % Smooth amplitude envelope using moving average
    %
    % Args:
    %   envelope: amplitude values
    %   window_size: number of samples to average
    %
    % Returns:
    %   smoothed: smoothed envelope

    if nargin < 2
        window_size = 5;
    end

    % Ensure envelope is a column vector
    envelope = envelope(:);

    % Ensure we have enough data
    if length(envelope) < window_size
        smoothed = envelope;
        return;
    end

    % Simple moving average using convolution
    kernel = ones(window_size, 1) / window_size;  % Make kernel same orientation
    smoothed = conv(envelope, kernel, 'same');
end


function [peak_times, peak_values] = find_amplitude_peaks(times, amps, min_peak_height, min_peak_distance)
    % Find peaks in amplitude trajectory
    %
    % Args:
    %   times: time vector (seconds)
    %   amps: amplitude values (dBFS)
    %   min_peak_height: minimum amplitude for peak (dBFS)
    %   min_peak_distance: minimum samples between peaks
    %
    % Returns:
    %   peak_times: time points of peaks
    %   peak_values: amplitude values at peaks

    if nargin < 3
        min_peak_height = -40;
    end
    if nargin < 4
        min_peak_distance = 10;
    end

    % Only consider valid (non-zero) amplitudes
    valid = amps > -60;

    if sum(valid) < 2
        peak_times = [];
        peak_values = [];
        return;
    end

    valid_times = times(valid);
    valid_amps = amps(valid);

    % Shift amplitudes to be positive for findpeaks
    min_amp = min(valid_amps);
    shifted_amps = valid_amps - min_amp;
    shifted_threshold = min_peak_height - min_amp;

    % Make sure threshold is non-negative
    if shifted_threshold < 0
        shifted_threshold = 0;
    end

    % Find peaks
    [peaks, locs] = findpeaks(shifted_amps, 'MinPeakHeight', shifted_threshold, ...
                                            'MinPeakDistance', min_peak_distance);

    if isempty(locs)
        peak_times = [];
        peak_values = [];
        return;
    end

    % Convert back to original dBFS scale
    peak_times = valid_times(locs);
    peak_values = peaks + min_amp;
end


function [rate_of_change, inflection_times] = analyze_envelope_dynamics(times, amps)
    % Analyze rate of change in amplitude envelope
    %
    % Args:
    %   times: time vector (seconds)
    %   amps: amplitude values (dBFS)
    %
    % Returns:
    %   rate_of_change: dB/second change rate
    %   inflection_times: times where derivative changes sign

    if length(times) < 2
        rate_of_change = [];
        inflection_times = [];
        return;
    end

    % First derivative (rate of change)
    dt = diff(times);
    damp = diff(amps);
    rate_of_change = damp ./ dt;

    % Smooth the rate of change to reduce noise
    if length(rate_of_change) > 10
        rate_of_change = smooth_envelope(rate_of_change, 10);
    end

    if length(rate_of_change) < 2
        inflection_times = [];
        return;
    end

    % Find inflection points (where derivative changes sign significantly)
    % Use a threshold to avoid detecting tiny fluctuations

    % Calculate a dynamic threshold based on the data
    threshold = max(abs(rate_of_change)) * 0.1;  % 10% of max rate
    if threshold < 5
        threshold = 5;  % Minimum threshold of 5 dB/s
    end

    % Only consider significant changes
    sign_roc = sign(rate_of_change);
    significant = abs(rate_of_change) > threshold;
    sign_roc(~significant) = 0;

    % Find sign changes
    sign_changes = diff(sign_roc);
    inflection_indices = find(abs(sign_changes) >= 2);  % Sign change from -1 to +1 or vice versa

    if ~isempty(inflection_indices)
        % Map back to time indices
        time_indices = inflection_indices + 1;
        time_indices = time_indices(time_indices <= length(times));

        % Further filter: ensure inflections are spaced apart
        if length(time_indices) > 1
            min_spacing = 0.1;  % At least 0.1 seconds between inflections
            filtered_indices = [time_indices(1)];
            for i = 2:length(time_indices)
                if times(time_indices(i)) - times(filtered_indices(end)) > min_spacing
                    filtered_indices = [filtered_indices; time_indices(i)];
                end
            end
            time_indices = filtered_indices;
        end

        inflection_times = times(time_indices);
    else
        inflection_times = [];
    end
end


function [period, confidence] = find_amplitude_periodicity(amps, fs_analysis)
    % Detect periodic patterns in amplitude using autocorrelation
    %
    % Args:
    %   amps: amplitude values (dBFS)
    %   fs_analysis: effective sample rate of amplitude data (Hz)
    %
    % Returns:
    %   period: detected period in seconds (0 if none)
    %   confidence: confidence in detection (0-1)

    % Only use valid samples
    valid = amps > -60;
    signal = amps(valid);

    if length(signal) < 10
        period = 0;
        confidence = 0;
        return;
    end

    % Autocorrelation
    [acf, lags] = xcorr(signal, 'coeff');

    % Look for peaks in positive lags only
    half = floor(length(acf)/2) + 1;
    acf_pos = acf(half:end);
    lags_pos = lags(half:end);

    % Shift ACF to be non-negative for findpeaks
    min_acf = min(acf_pos);
    if min_acf < 0
        acf_shifted = acf_pos - min_acf;
    else
        acf_shifted = acf_pos;
    end

    % Find first significant peak after lag 0
    peak_threshold = 0.3 * (max(acf_shifted) - min(acf_shifted));

    [peaks, locs] = findpeaks(acf_shifted, 'MinPeakHeight', peak_threshold, ...
                                           'MinPeakDistance', 5);

    if ~isempty(locs) && locs(1) > 1
        period = lags_pos(locs(1)) / fs_analysis;
        confidence = acf_pos(locs(1));
    else
        period = 0;
        confidence = 0;
    end
end


function segments = segment_by_activity(times, amps, threshold, min_duration)
    % Segment envelope into active and quiet regions
    %
    % Args:
    %   times: time vector (seconds)
    %   amps: amplitude values (dBFS)
    %   threshold: amplitude threshold for "active" (dBFS)
    %   min_duration: minimum segment duration (seconds)
    %
    % Returns:
    %   segments: array of structs with start, end_time, duration, mean_amp, max_amp

    if nargin < 3
        threshold = -40;
    end
    if nargin < 4
        min_duration = 0.1;
    end

    % Find active regions
    active = amps > threshold;

    % Find transitions
    transitions = diff([0; active; 0]);
    starts = find(transitions == 1);
    ends = find(transitions == -1) - 1;

    segments = [];
    for i = 1:length(starts)
        if starts(i) <= length(times) && ends(i) <= length(times)
            duration = times(ends(i)) - times(starts(i));
            if duration >= min_duration
                segment.start = times(starts(i));
                segment.end_time = times(ends(i));
                segment.duration = duration;
                segment.mean_amp = mean(amps(starts(i):ends(i)));
                segment.max_amp = max(amps(starts(i):ends(i)));
                segments = [segments; segment];
            end
        end
    end
end


function analysis = analyze_partial_rhythm(times, amps, partial_idx, options)
    % Comprehensive rhythmic analysis of a partial's amplitude envelope
    %
    % Args:
    %   times: time vector (seconds)
    %   amps: amplitude trajectory for this partial (dBFS)
    %   partial_idx: which partial this is (for labeling)
    %   options: struct with analysis parameters (optional)
    %     .smooth_window - samples for smoothing (default 5)
    %     .activity_threshold - dBFS threshold for activity (default -40)
    %     .min_segment_duration - minimum segment length in seconds (default 0.1)
    %
    % Returns:
    %   analysis: struct with rhythmic features

    if nargin < 4
        options = struct();
    end

    % Default options
    if ~isfield(options, 'smooth_window')
        options.smooth_window = 5;
    end
    if ~isfield(options, 'activity_threshold')
        options.activity_threshold = -40;
    end
    if ~isfield(options, 'min_segment_duration')
        options.min_segment_duration = 0.1;
    end

    analysis.partial_idx = partial_idx;

    % 1. Basic statistics
    valid = amps > -60;

    if sum(valid) == 0
        % No valid data for this partial
        analysis.mean_amplitude = -60;
        analysis.max_amplitude = -60;
        analysis.min_amplitude = -60;
        analysis.amplitude_range = 0;
        analysis.smoothed = amps;
        analysis.num_peaks = 0;
        analysis.peak_times = [];
        analysis.peak_amplitudes = [];
        analysis.inter_peak_intervals = [];
        analysis.num_segments = 0;
        analysis.segments = [];
        analysis.rate_of_change = [];
        analysis.inflection_points = [];
        analysis.num_inflections = 0;
        analysis.mean_rate_of_change = 0;
        analysis.detected_period = 0;
        analysis.periodicity_confidence = 0;
        analysis.envelope_type = 'silent';
        analysis.mean_rhythm = 0;
        analysis.rhythm_regularity = 0;
        analysis.min_interval = 0;
        analysis.max_interval = 0;
        analysis.total_duration = times(end) - times(1);
        analysis.mean_segment_duration = 0;
        analysis.total_active_time = 0;
        analysis.activity_ratio = 0;
        return;
    end

    analysis.mean_amplitude = mean(amps(valid));
    analysis.max_amplitude = max(amps(valid));
    analysis.min_amplitude = min(amps(valid));
    analysis.amplitude_range = analysis.max_amplitude - analysis.min_amplitude;

    % 2. Smoothed envelope
    analysis.smoothed = smooth_envelope(amps, options.smooth_window);

    % 3. Find peaks on smoothed data (reduces false peaks)
    [peak_times, peak_amps] = find_amplitude_peaks(times, analysis.smoothed, ...
                                                    analysis.mean_amplitude, ...
                                                    20);  % Increased min distance
    analysis.peak_times = peak_times;
    analysis.peak_amplitudes = peak_amps;
    analysis.num_peaks = length(peak_times);

    % 4. Inter-peak intervals (rhythmic spacing)
    if length(peak_times) > 1
        analysis.inter_peak_intervals = diff(peak_times);
        analysis.mean_rhythm = mean(analysis.inter_peak_intervals);
        analysis.rhythm_regularity = std(analysis.inter_peak_intervals);
        analysis.min_interval = min(analysis.inter_peak_intervals);
        analysis.max_interval = max(analysis.inter_peak_intervals);
    else
        analysis.inter_peak_intervals = [];
        analysis.mean_rhythm = 0;
        analysis.rhythm_regularity = 0;
        analysis.min_interval = 0;
        analysis.max_interval = 0;
    end

    % 5. Activity segments
    analysis.segments = segment_by_activity(times, amps, ...
                                            options.activity_threshold, ...
                                            options.min_segment_duration);
    analysis.num_segments = length(analysis.segments);

    % 6. Dynamics (rate of change) - use smoothed data
    valid_times = times(valid);
    valid_amps = analysis.smoothed(valid);

    if length(valid_times) > 2
        [roc, inflections] = analyze_envelope_dynamics(valid_times, valid_amps);
        analysis.rate_of_change = roc;
        analysis.inflection_points = inflections;
        analysis.num_inflections = length(inflections);
        analysis.mean_rate_of_change = mean(abs(roc));
    else
        analysis.rate_of_change = [];
        analysis.inflection_points = [];
        analysis.num_inflections = 0;
        analysis.mean_rate_of_change = 0;
    end

    % 7. Periodicity analysis
    if sum(valid) > 10
        fs_analysis = 1 / mean(diff(times));
        [period, conf] = find_amplitude_periodicity(amps, fs_analysis);
        analysis.detected_period = period;
        analysis.periodicity_confidence = conf;
    else
        analysis.detected_period = 0;
        analysis.periodicity_confidence = 0;
    end

    % 8. Envelope shape classification
    if analysis.num_segments > 0
        first_seg = analysis.segments(1);
        if first_seg.mean_amp > analysis.mean_amplitude
            analysis.envelope_type = 'percussive';
        else
            analysis.envelope_type = 'sustained';
        end
    else
        analysis.envelope_type = 'sparse';
    end

    % 9. Duration statistics
    analysis.total_duration = times(end) - times(1);
    if analysis.num_segments > 0
        segment_durations = arrayfun(@(s) s.duration, analysis.segments);
        analysis.mean_segment_duration = mean(segment_durations);
        analysis.total_active_time = sum(segment_durations);
        analysis.activity_ratio = analysis.total_active_time / analysis.total_duration;
    else
        analysis.mean_segment_duration = 0;
        analysis.total_active_time = 0;
        analysis.activity_ratio = 0;
    end
end

% ============================================================================
% RHYTHMIC ANALYSIS VISUALIZATION FUNCTIONS
% ============================================================================

function plot_rhythmic_analysis(times, amps, analysis, title_str, time_range)
    % Visualize rhythmic analysis of a single partial
    %
    % Args:
    %   times: time vector (seconds)
    %   amps: amplitude trajectory (dBFS)
    %   analysis: analysis struct from analyze_partial_rhythm
    %   title_str: optional title (default '')
    %   time_range: [min, max] time range to plot (default: auto)

    if nargin < 4
        title_str = '';
    end
    if nargin < 5
        time_range = [];
    end

    figure;

    % Subplot 1: Amplitude envelope with peaks and segments
    subplot(3, 1, 1);
    hold on;

    % Plot raw envelope
    valid = amps > -60;
    plot(times(valid), amps(valid), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);

    % Plot smoothed envelope
    plot(times(valid), analysis.smoothed(valid), '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2);

    % Mark peaks
    if ~isempty(analysis.peak_times)
        plot(analysis.peak_times, analysis.peak_amplitudes, 'ro', ...
             'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end

    % Shade active segments
    if analysis.num_segments > 0
        ylims = ylim;
        for i = 1:length(analysis.segments)
            seg = analysis.segments(i);
            fill([seg.start seg.end_time seg.end_time seg.start], ...
                 [ylims(1) ylims(1) ylims(2) ylims(2)], ...
                 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
    end

    hold off;

    if ~isempty(time_range)
        xlim(time_range);
    end

    ylabel('Amplitude (dBFS)');
    if ~isempty(title_str)
        title(sprintf('%s - Partial %d Envelope', title_str, analysis.partial_idx));
    else
        title(sprintf('Partial %d Envelope', analysis.partial_idx));
    end
    legend({'Raw', 'Smoothed', 'Peaks', 'Active'}, 'location', 'eastoutside');
    grid on;

    % Subplot 2: Rate of change (dynamics)
    subplot(3, 1, 2);
    if ~isempty(analysis.rate_of_change)
        valid_times = times(valid);
        plot(valid_times(1:length(analysis.rate_of_change)), analysis.rate_of_change, '-', ...
             'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);
        hold on;

        % Mark inflection points
        if ~isempty(analysis.inflection_points)
            % Get y-values at inflection points
            y_vals = zeros(size(analysis.inflection_points));
            for i = 1:length(analysis.inflection_points)
                [~, idx] = min(abs(valid_times(1:length(analysis.rate_of_change)) - analysis.inflection_points(i)));
                if idx <= length(analysis.rate_of_change)
                    y_vals(i) = analysis.rate_of_change(idx);
                end
            end
            plot(analysis.inflection_points, y_vals, ...
                 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
        end

        % Zero line
        plot(xlim, [0 0], 'k--', 'LineWidth', 0.5);
        hold off;

        if ~isempty(time_range)
            xlim(time_range);
        end

        ylabel('Rate of Change (dB/s)');
        title(sprintf('Amplitude Dynamics (%d inflections)', analysis.num_inflections));
        if ~isempty(analysis.inflection_points)
            legend({'dA/dt', 'Inflections', 'Zero'}, 'location', 'eastoutside');
        else
            legend({'dA/dt', 'Zero'}, 'location', 'eastoutside');
        end
        grid on;
    else
        text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center');
        axis off;
    end

    % Subplot 3: Inter-peak intervals (rhythm)
    subplot(3, 1, 3);
    if ~isempty(analysis.inter_peak_intervals) && length(analysis.inter_peak_intervals) > 0
        bar(1:length(analysis.inter_peak_intervals), analysis.inter_peak_intervals, ...
            'FaceColor', [0.4 0.6 0.4]);
        hold on;

        % Mean line
        plot(xlim, [analysis.mean_rhythm analysis.mean_rhythm], 'r--', ...
             'LineWidth', 2);
        hold off;

        xlabel('Peak Pair Index');
        ylabel('Interval (seconds)');
        title(sprintf('Inter-Peak Intervals (Mean: %.3f s, Std: %.3f s)', ...
                     analysis.mean_rhythm, analysis.rhythm_regularity));
        legend({'Intervals', 'Mean'}, 'location', 'northeast');
        grid on;
    else
        text(0.5, 0.5, sprintf('Insufficient peaks (%d found)', analysis.num_peaks), ...
             'HorizontalAlignment', 'center');
        axis off;
    end
end


function plot_rhythmic_summary(analyses, title_str)
    % Create summary visualization comparing all partials
    %
    % Args:
    %   analyses: cell array of analysis structs
    %   title_str: optional title (default '')

    if nargin < 2
        title_str = '';
    end

    num_partials = length(analyses);

    figure;

    % Subplot 1: Mean rhythm comparison
    subplot(2, 2, 1);
    mean_rhythms = zeros(num_partials, 1);
    for i = 1:num_partials
        mean_rhythms(i) = analyses{i}.mean_rhythm;
    end
    bar(1:num_partials, mean_rhythms, 'FaceColor', [0.3 0.5 0.7]);
    xlabel('Partial Number');
    ylabel('Mean Rhythm (seconds)');
    title('Mean Inter-Peak Interval by Partial');
    grid on;

    % Subplot 2: Number of peaks
    subplot(2, 2, 2);
    num_peaks = zeros(num_partials, 1);
    for i = 1:num_partials
        num_peaks(i) = analyses{i}.num_peaks;
    end
    bar(1:num_partials, num_peaks, 'FaceColor', [0.7 0.3 0.3]);
    xlabel('Partial Number');
    ylabel('Number of Peaks');
    title('Rhythmic Density by Partial');
    grid on;

    % Subplot 3: Activity ratio
    subplot(2, 2, 3);
    activity_ratios = zeros(num_partials, 1);
    for i = 1:num_partials
        activity_ratios(i) = analyses{i}.activity_ratio * 100;
    end
    bar(1:num_partials, activity_ratios, 'FaceColor', [0.3 0.7 0.3]);
    xlabel('Partial Number');
    ylabel('Activity (%)');
    title('Time Active by Partial');
    ylim([0 100]);
    grid on;

    % Subplot 4: Periodicity confidence
    subplot(2, 2, 4);
    periodicities = zeros(num_partials, 1);
    confidences = zeros(num_partials, 1);
    for i = 1:num_partials
        periodicities(i) = analyses{i}.detected_period;
        confidences(i) = analyses{i}.periodicity_confidence;
    end

    % Only plot partials with detected periods
    has_period = periodicities > 0;
    if sum(has_period) > 0
        scatter(find(has_period), periodicities(has_period), ...
                100 * confidences(has_period) + 20, ...
                confidences(has_period), ...
                'filled');
        colorbar;
        xlabel('Partial Number');
        ylabel('Detected Period (seconds)');
        title('Periodicity Detection (size = confidence)');
        grid on;
    else
        text(0.5, 0.5, 'No periodicities detected', ...
             'HorizontalAlignment', 'center');
        axis off;
    end

    % Overall title
    if ~isempty(title_str)
        annotation('textbox', [0 0.96 1 0.04], ...
                   'String', [title_str ' - Rhythmic Analysis Summary'], ...
                   'EdgeColor', 'none', ...
                   'HorizontalAlignment', 'center', ...
                   'FontSize', 12, ...
                   'FontWeight', 'bold');
    end
end


function plot_rhythm_matrix(analyses, title_str)
    % Create a heat map showing rhythmic patterns across all partials over time
    %
    % Args:
    %   analyses: cell array of analysis structs
    %   title_str: optional title (default '')

    if nargin < 2
        title_str = '';
    end

    num_partials = length(analyses);

    figure;

    % Find global time range
    max_time = 0;
    for i = 1:num_partials
        if ~isempty(analyses{i}.peak_times)
            max_time = max(max_time, max(analyses{i}.peak_times));
        end
    end

    if max_time == 0
        text(0.5, 0.5, 'No rhythmic data to display', ...
             'HorizontalAlignment', 'center');
        return;
    end

    % Create time bins (100 bins across duration)
    num_bins = 100;
    time_bins = linspace(0, max_time, num_bins);
    bin_width = time_bins(2) - time_bins(1);

    % Create activity matrix
    activity_matrix = zeros(num_partials, num_bins);

    for p = 1:num_partials
        if ~isempty(analyses{p}.peak_times)
            % Bin the peaks
            for peak_time = analyses{p}.peak_times'
                bin_idx = min(floor(peak_time / bin_width) + 1, num_bins);
                activity_matrix(p, bin_idx) = activity_matrix(p, bin_idx) + 1;
            end
        end
    end

    % Plot heat map
    imagesc(time_bins, 1:num_partials, activity_matrix);
    colormap(hot);
    colorbar;

    xlabel('Time (seconds)');
    ylabel('Partial Number');
    if ~isempty(title_str)
        title([title_str ' - Rhythmic Activity Heat Map']);
    else
        title('Rhythmic Activity Heat Map');
    end

    % Set y-axis to show all partials
    yticks(1:num_partials);
    set(gca, 'YDir', 'normal');
end

printf('Rhythmic analysis functions loaded!\n');
