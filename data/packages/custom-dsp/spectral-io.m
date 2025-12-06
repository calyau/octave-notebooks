% spectral-io.m
% Functions for saving and loading spectral analysis data in .mat format

function save_static_spectrum(filename, freqs, amps, metadata)
    % Save a single spectral snapshot to .mat file
    %
    % Parameters:
    %   filename - output .mat file path (string)
    %   freqs - 1D array of frequencies in Hz
    %   amps - 1D array of amplitudes in dBFS
    %   metadata - struct with fields:
    %     .fundamental (required, scalar) - fundamental frequency in Hz
    %     .sample_rate (required, scalar) - sample rate in Hz
    %     .instrument_id (optional, scalar, default 0) - numeric instrument code
    %     .midi_note (optional, scalar, default 0) - MIDI note number (0-127)
    %     .time_stamp (optional, scalar, default 0.0) - time in seconds

    % Reshape to row vectors
    freqs = reshape(freqs, 1, []);
    amps = reshape(amps, 1, []);

    % Calculate number of partials
    num_partials = length(freqs);

    % Extract required metadata
    fundamental = metadata.fundamental;
    sample_rate = metadata.sample_rate;

    % Extract optional metadata with defaults
    if isfield(metadata, 'instrument_id')
        instrument_id = metadata.instrument_id;
    else
        instrument_id = 0;
    end

    if isfield(metadata, 'midi_note')
        midi_note = metadata.midi_note;
    else
        midi_note = 0;
    end

    if isfield(metadata, 'time_stamp')
        time_stamp = metadata.time_stamp;
    else
        time_stamp = 0.0;
    end

    % Format signature for mat2sdif tool detection
    octave_dsp_format = 'spectral_static';
    octave_dsp_version = 1.0;

    % Save to .mat file (MATLAB v5/v6 compatible)
    save('-v6', filename, 'freqs', 'amps', 'num_partials', ...
         'fundamental', 'sample_rate', 'time_stamp', ...
         'instrument_id', 'midi_note', ...
         'octave_dsp_format', 'octave_dsp_version');

    % Print confirmation
    printf('Saved static spectrum: %s (%d partials)\n', filename, num_partials);
end


function save_spectral_evolution(filename, times, freqs_matrix, amps_matrix, metadata)
    % Save time-varying spectral data to .mat file
    %
    % Parameters:
    %   filename - output .mat file path (string)
    %   times - 1D array of time points in seconds
    %   freqs_matrix - 2D array [num_frames × max_partials], zero-padded
    %   amps_matrix - 2D array [num_frames × max_partials], zero-padded
    %   metadata - struct with fields:
    %     .num_partials_per_frame (required, 1D array) - actual number of partials in each frame
    %     .fundamental (required, scalar) - fundamental frequency in Hz
    %     .sample_rate (required, scalar) - sample rate in Hz
    %     .window_size (required, scalar) - FFT window size in samples
    %     .hop_size (required, scalar) - hop size in samples between frames
    %     .instrument_id (optional, scalar, default 0)
    %     .midi_note (optional, scalar, default 0)

    % Reshape times to row vector
    times = reshape(times, 1, []);

    % Calculate number of frames
    num_frames = length(times);

    % Extract required metadata
    num_partials_per_frame = reshape(metadata.num_partials_per_frame, 1, []);
    fundamental = metadata.fundamental;
    sample_rate = metadata.sample_rate;
    window_size = metadata.window_size;
    hop_size = metadata.hop_size;

    % Extract optional metadata with defaults
    if isfield(metadata, 'instrument_id')
        instrument_id = metadata.instrument_id;
    else
        instrument_id = 0;
    end

    if isfield(metadata, 'midi_note')
        midi_note = metadata.midi_note;
    else
        midi_note = 0;
    end

    % Rename matrices for saving
    freqs = freqs_matrix;
    amps = amps_matrix;

    % Format signature for mat2sdif tool detection
    octave_dsp_format = 'spectral_evolution';
    octave_dsp_version = 1.0;

    % Save to .mat file (MATLAB v5/v6 compatible)
    save('-v6', filename, 'times', 'freqs', 'amps', ...
         'num_frames', 'num_partials_per_frame', ...
         'fundamental', 'sample_rate', 'window_size', 'hop_size', ...
         'instrument_id', 'midi_note', ...
         'octave_dsp_format', 'octave_dsp_version');

    % Print confirmation
    printf('Saved spectral evolution: %s (%d frames, %.2f-%.2f sec)\n', ...
           filename, num_frames, times(1), times(end));
end

printf('Spectral I/O functions loaded!\n');
