% ============================================================================
% RHYTHMIC ANALYSIS EXPORT FUNCTIONS
% ============================================================================

function export_rhythmic_analysis_text(analyses, filename, title_str)
    % Export rhythmic analysis to human-readable text file
    %
    % Args:
    %   analyses: cell array of analysis structs
    %   filename: output text file path
    %   title_str: optional title for the report

    if nargin < 3
        title_str = 'Rhythmic Analysis Report';
    end

    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file for writing: %s', filename);
    end

    % Write header
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '%s\n', title_str);
    fprintf(fid, 'Generated: %s\n', datestr(now));
    fprintf(fid, '================================================================================\n\n');

    % Write summary statistics
    fprintf(fid, 'OVERALL SUMMARY\n');
    fprintf(fid, '---------------\n');
    fprintf(fid, 'Number of partials analyzed: %d\n\n', length(analyses));

    % Write individual partial analyses
    for i = 1:length(analyses)
        a = analyses{i};

        fprintf(fid, '================================================================================\n');
        fprintf(fid, 'PARTIAL %d\n', a.partial_idx);
        fprintf(fid, '================================================================================\n\n');

        fprintf(fid, 'Amplitude Characteristics:\n');
        fprintf(fid, '  Mean amplitude:     %6.1f dBFS\n', a.mean_amplitude);
        fprintf(fid, '  Max amplitude:      %6.1f dBFS\n', a.max_amplitude);
        fprintf(fid, '  Min amplitude:      %6.1f dBFS\n', a.min_amplitude);
        fprintf(fid, '  Amplitude range:    %6.1f dB\n', a.amplitude_range);
        fprintf(fid, '  Envelope type:      %s\n\n', a.envelope_type);

        fprintf(fid, 'Rhythmic Features:\n');
        fprintf(fid, '  Number of peaks:    %d\n', a.num_peaks);
        if a.num_peaks > 1
            fprintf(fid, '  Mean rhythm:        %.4f seconds (%.2f Hz)\n', ...
                   a.mean_rhythm, 1/a.mean_rhythm);
            fprintf(fid, '  Rhythm regularity:  %.4f seconds (std dev)\n', a.rhythm_regularity);
            fprintf(fid, '  Min interval:       %.4f seconds\n', a.min_interval);
            fprintf(fid, '  Max interval:       %.4f seconds\n', a.max_interval);

            % Write all intervals
            fprintf(fid, '\n  Inter-peak intervals:\n');
            for j = 1:length(a.inter_peak_intervals)
                fprintf(fid, '    [%2d] %.4f seconds\n', j, a.inter_peak_intervals(j));
            end
        else
            fprintf(fid, '  (insufficient peaks for rhythm analysis)\n');
        end
        fprintf(fid, '\n');

        fprintf(fid, 'Activity Segmentation:\n');
        fprintf(fid, '  Number of segments: %d\n', a.num_segments);
        fprintf(fid, '  Total duration:     %.3f seconds\n', a.total_duration);
        fprintf(fid, '  Active time:        %.3f seconds (%.1f%%)\n', ...
               a.total_active_time, a.activity_ratio * 100);
        if a.num_segments > 0
            fprintf(fid, '  Mean segment:       %.3f seconds\n\n', a.mean_segment_duration);
            fprintf(fid, '  Segment details:\n');
            for j = 1:length(a.segments)
                seg = a.segments(j);
                fprintf(fid, '    [%2d] %.3f - %.3f s (dur: %.3f s, mean: %.1f dBFS, max: %.1f dBFS)\n', ...
                       j, seg.start, seg.end_time, seg.duration, seg.mean_amp, seg.max_amp);
            end
        end
        fprintf(fid, '\n');

        fprintf(fid, 'Dynamics:\n');
        fprintf(fid, '  Inflection points:  %d\n', a.num_inflections);
        fprintf(fid, '  Mean rate of change: %.2f dB/s\n\n', a.mean_rate_of_change);

        fprintf(fid, 'Periodicity:\n');
        if a.detected_period > 0
            fprintf(fid, '  Detected period:    %.4f seconds (%.2f Hz)\n', ...
                   a.detected_period, 1/a.detected_period);
            fprintf(fid, '  Confidence:         %.2f\n', a.periodicity_confidence);
        else
            fprintf(fid, '  No periodicity detected\n');
        end

        fprintf(fid, '\n\n');
    end

    fclose(fid);
    printf('Rhythmic analysis exported to: %s\n', filename);
end


function export_rhythmic_analysis_csv(analyses, filename)
    % Export rhythmic analysis to CSV for spreadsheet/data analysis
    %
    % Args:
    %   analyses: cell array of analysis structs
    %   filename: output CSV file path

    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file for writing: %s', filename);
    end

    % Write header
    fprintf(fid, 'partial,envelope_type,mean_amp,max_amp,min_amp,amp_range,');
    fprintf(fid, 'num_peaks,mean_rhythm,rhythm_regularity,min_interval,max_interval,');
    fprintf(fid, 'num_segments,total_duration,active_time,activity_ratio,mean_segment_duration,');
    fprintf(fid, 'num_inflections,mean_rate_of_change,detected_period,periodicity_confidence\n');

    % Write data rows
    for i = 1:length(analyses)
        a = analyses{i};

        fprintf(fid, '%d,%s,%.2f,%.2f,%.2f,%.2f,', ...
               a.partial_idx, a.envelope_type, ...
               a.mean_amplitude, a.max_amplitude, a.min_amplitude, a.amplitude_range);

        fprintf(fid, '%d,%.4f,%.4f,%.4f,%.4f,', ...
               a.num_peaks, a.mean_rhythm, a.rhythm_regularity, ...
               a.min_interval, a.max_interval);

        fprintf(fid, '%d,%.4f,%.4f,%.4f,%.4f,', ...
               a.num_segments, a.total_duration, a.total_active_time, ...
               a.activity_ratio, a.mean_segment_duration);

        fprintf(fid, '%d,%.4f,%.4f,%.4f\n', ...
               a.num_inflections, a.mean_rate_of_change, ...
               a.detected_period, a.periodicity_confidence);
    end

    fclose(fid);
    printf('Rhythmic analysis exported to CSV: %s\n', filename);
end


function export_rhythm_for_notation(analyses, filename, tempo_bpm)
    % Export rhythmic patterns in a format useful for notation software
    % Converts time intervals to musical durations at given tempo
    %
    % Args:
    %   analyses: cell array of analysis structs
    %   filename: output text file path
    %   tempo_bpm: tempo in beats per minute (default 60)

    if nargin < 3
        tempo_bpm = 60;
    end

    % Calculate beat duration in seconds
    beat_duration = 60 / tempo_bpm;

    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file for writing: %s', filename);
    end

    fprintf(fid, 'RHYTHMIC NOTATION GUIDE\n');
    fprintf(fid, '=======================\n');
    fprintf(fid, 'Tempo: %d BPM (quarter note = %.3f seconds)\n\n', tempo_bpm, beat_duration);

    for i = 1:length(analyses)
        a = analyses{i};

        if a.num_peaks < 2
            continue;  % Skip partials without rhythm
        end

        fprintf(fid, '\n--- PARTIAL %d ---\n', a.partial_idx);
        fprintf(fid, 'Mean interval: %.4f seconds = %.2f beats\n\n', ...
               a.mean_rhythm, a.mean_rhythm / beat_duration);

        fprintf(fid, 'Rhythmic sequence (in beats):\n');

        for j = 1:length(a.inter_peak_intervals)
            interval_sec = a.inter_peak_intervals(j);
            interval_beats = interval_sec / beat_duration;

            % Suggest closest standard note values
            note_suggestion = suggest_note_value(interval_beats);

            fprintf(fid, '  %2d: %.4f s = %.3f beats ≈ %s\n', ...
                   j, interval_sec, interval_beats, note_suggestion);
        end

        % Suggest a rhythmic pattern
        fprintf(fid, '\nSuggested pattern: ');
        for j = 1:min(8, length(a.inter_peak_intervals))  % First 8 intervals
            interval_beats = a.inter_peak_intervals(j) / beat_duration;
            fprintf(fid, '%s ', suggest_note_value(interval_beats));
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
    printf('Notation guide exported to: %s\n', filename);
end


function note_str = suggest_note_value(beats)
    % Helper function to suggest closest standard note value
    %
    % Args:
    %   beats: duration in quarter note beats
    %
    % Returns:
    %   note_str: string describing the note value

    % Standard note values in beats
    note_values = [4, 3, 2, 1.5, 1, 0.75, 0.5, 0.375, 0.25, 0.1875, 0.125];
    note_names = {'whole', 'dotted half', 'half', 'dotted quarter', ...
                  'quarter', 'dotted eighth', 'eighth', 'dotted sixteenth', ...
                  'sixteenth', 'dotted 32nd', '32nd'};

    % Find closest match
    [~, idx] = min(abs(note_values - beats));
    note_str = note_names{idx};

    % Add deviation if significant
    deviation = beats - note_values(idx);
    if abs(deviation) > 0.1
        if deviation > 0
            note_str = [note_str ' (+)'];
        else
            note_str = [note_str ' (-)'];
        end
    end
end


function export_rhythm_json(analyses, filename, metadata)
    % Export rhythmic analysis to JSON format for use in other tools
    %
    % Args:
    %   analyses: cell array of analysis structs
    %   filename: output JSON file path
    %   metadata: optional struct with additional metadata

    if nargin < 3
        metadata = struct();
    end

    % Build JSON structure
    output = struct();
    output.metadata = metadata;
    output.metadata.export_date = datestr(now);
    output.metadata.num_partials = length(analyses);
    output.partials = analyses;

    % Convert to JSON string (manual since Octave might not have jsonencode)
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file for writing: %s', filename);
    end

    fprintf(fid, '{\n');
    fprintf(fid, '  "metadata": {\n');
    fprintf(fid, '    "export_date": "%s",\n', datestr(now));
    fprintf(fid, '    "num_partials": %d\n', length(analyses));
    fprintf(fid, '  },\n');
    fprintf(fid, '  "partials": [\n');

    for i = 1:length(analyses)
        a = analyses{i};

        fprintf(fid, '    {\n');
        fprintf(fid, '      "partial_idx": %d,\n', a.partial_idx);
        fprintf(fid, '      "envelope_type": "%s",\n', a.envelope_type);
        fprintf(fid, '      "mean_amplitude": %.2f,\n', a.mean_amplitude);
        fprintf(fid, '      "amplitude_range": %.2f,\n', a.amplitude_range);
        fprintf(fid, '      "num_peaks": %d,\n', a.num_peaks);
        fprintf(fid, '      "mean_rhythm": %.4f,\n', a.mean_rhythm);
        fprintf(fid, '      "rhythm_regularity": %.4f,\n', a.rhythm_regularity);
        fprintf(fid, '      "activity_ratio": %.4f,\n', a.activity_ratio);
        fprintf(fid, '      "detected_period": %.4f,\n', a.detected_period);
        fprintf(fid, '      "periodicity_confidence": %.4f,\n', a.periodicity_confidence);

        % Export peak times
        fprintf(fid, '      "peak_times": [');
        if ~isempty(a.peak_times)
            for j = 1:length(a.peak_times)
                fprintf(fid, '%.4f', a.peak_times(j));
                if j < length(a.peak_times)
                    fprintf(fid, ', ');
                end
            end
        end
        fprintf(fid, '],\n');

        % Export inter-peak intervals
        fprintf(fid, '      "inter_peak_intervals": [');
        if ~isempty(a.inter_peak_intervals)
            for j = 1:length(a.inter_peak_intervals)
                fprintf(fid, '%.4f', a.inter_peak_intervals(j));
                if j < length(a.inter_peak_intervals)
                    fprintf(fid, ', ');
                end
            end
        end
        fprintf(fid, ']\n');

        fprintf(fid, '    }');
        if i < length(analyses)
            fprintf(fid, ',');
        end
        fprintf(fid, '\n');
    end

    fprintf(fid, '  ]\n');
    fprintf(fid, '}\n');

    fclose(fid);
    printf('Rhythmic analysis exported to JSON: %s\n', filename);
end

printf('Rhythmic I/O functions loaded ✓\n');
