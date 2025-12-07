% Find harmonics from a spectrum
function [freqs, amps, metadata] = analyze_static_spectrum(signal, fs, window_size)
    % Analyze a single spectral snapshot
    %
    % Args:
    %   signal: audio signal (mono)
    %   fs: sample rate in Hz
    %   window_size: FFT window size (default 4096)
    %
    % Returns:
    %   freqs: array of peak frequencies in Hz
    %   amps: array of peak amplitudes in dBFS
    %   metadata: struct with analysis parameters

    if nargin < 3
        window_size = 4096;
    end

    % Take a chunk from the middle of the file (avoid silence at start/end)
    mid_point = floor(length(signal) / 2);
    start_idx = mid_point - floor(window_size / 2);
    if start_idx < 1
        start_idx = 1;
    end
    chunk_size = min(window_size, length(signal) - start_idx + 1);
    chunk = signal(start_idx:start_idx + chunk_size - 1) .* hanning(chunk_size);

    % Zero-pad for better frequency resolution
    padded = [chunk; zeros(window_size - chunk_size, 1)];

    % Compute spectrum with proper dBFS scaling
    spectrum = abs(fft(padded));
    spectrum = spectrum / (window_size / 2);  % Normalize by FFT size
    spectrum_db = 20 * log10(spectrum(1:window_size/2) + 1e-10);

    % Shift spectrum to be positive for findpeaks
    min_db = min(spectrum_db);
    spectrum_shifted = spectrum_db - min_db;

    % Frequency axis
    freqs_full = (0:window_size/2-1)' * fs / window_size;

    % Calculate peak threshold (handle edge cases)
    peak_threshold = max(spectrum_shifted) - 42;
    if peak_threshold < 0
        peak_threshold = 0;
    end

    % Find peaks
    [pks, locs] = findpeaks(spectrum_shifted, 'MinPeakHeight', peak_threshold, ...
                                              'MinPeakDistance', 10);

    % Return actual dB values (not shifted)
    freqs = freqs_full(locs);
    amps = spectrum_db(locs);

    % Sort by frequency
    [freqs, idx] = sort(freqs);
    amps = amps(idx);

    % Build metadata struct
    metadata.num_partials = length(freqs);
    metadata.sample_rate = fs;
    metadata.window_size = window_size;
    metadata.time_stamp = 0.0;  % Single snapshot at t=0

    % Set fundamental automatically (first frequency)
    if ~isempty(freqs)
        metadata.fundamental = freqs(1);
        metadata.midi_note = freq_to_midi_note(metadata.fundamental);
    else
        metadata.fundamental = 0;  % No partials found
        metadata.midi_note = 'N/A';
    end
end

% Convert frequency to MIDI note name (e.g., "A4", "C#5")
function note_name = freq_to_midi_note(freq)
    if freq <= 0
        note_name = 'N/A';
        return;
    end

    % MIDI note 69 = A4 = 440 Hz
    % MIDI formula: n = 69 + 12*log2(f/440)
    midi_number = 69 + 12 * log2(freq / 440);
    midi_rounded = round(midi_number);

    % Note names
    note_names = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'};

    % Calculate octave (MIDI 0 = C-1, so octave = floor(midi/12) - 1)
    octave = floor(midi_rounded / 12) - 1;
    note_idx = mod(midi_rounded, 12) + 1;

    note_name = sprintf('%s%d', note_names{note_idx}, octave);
end

% Convert dBFS amplitude to MIDI velocity (0-127)
% Uses perceptual loudness mapping where dBFS range maps to velocity
function velocity = amplitude_to_midi(db_amplitude)
    % dBFS typically ranges from about -60 dB (very quiet) to 0 dB (full scale)
    % Map this to MIDI velocity 0-127 using a perceptual curve
    %
    % Key reference points based on typical dynamic range:
    % -60 dBFS or below -> velocity 1 (threshold of audibility)
    % -40 dBFS -> velocity ~20 (ppp)
    % -30 dBFS -> velocity ~40 (pp)
    % -20 dBFS -> velocity ~64 (mf)
    % -10 dBFS -> velocity ~96 (f)
    % -3 dBFS -> velocity 115 (ff)
    % 0 dBFS -> velocity 127 (fff)

    if db_amplitude <= -60
        velocity = 1;
    elseif db_amplitude >= 0
        velocity = 127;
    else
        % Use exponential mapping that approximates perception
        % This creates a curve where quieter sounds get more resolution
        normalized = (db_amplitude + 60) / 60;  % Map -60..0 to 0..1
        velocity = round(1 + 126 * (normalized ^ 0.7));  % Power law for perception
        velocity = max(1, min(127, velocity));  % Clamp to valid range
    end
end

% Convert MIDI velocity to dynamic marking (pppp to fff)
function dynamic = midi_to_dynamic(velocity)
    % Standard dynamic markings with MIDI velocity ranges
    % Based on typical performance practice and perception
    if velocity <= 0
        dynamic = '----';
    elseif velocity <= 16
        dynamic = 'pppp';
    elseif velocity <= 33
        dynamic = 'ppp';
    elseif velocity <= 49
        dynamic = 'pp';
    elseif velocity <= 64
        dynamic = 'p';
    elseif velocity <= 80
        dynamic = 'mp';
    elseif velocity <= 96
        dynamic = 'mf';
    elseif velocity <= 112
        dynamic = 'f';
    elseif velocity <= 120
        dynamic = 'ff';
    else
        dynamic = 'fff';
    end
end

% Display harmonics from a spectrum with enhanced MIDI information
function print_harmonics(peak_freqs, peak_amps)
    if isempty(peak_freqs)
        printf('No peaks found\n');
        return;
    end

    fundamental = peak_freqs(1);

    % Print header
    printf('%-12s %-10s %-10s %-13s %-12s %-10s %-8s\n', ...
           'Frequency', 'Closest', 'MIDI', 'Partial', 'Amplitude', 'Velocity', 'Dynamic');
    printf('%-12s %-10s %-10s %-13s %-12s %-10s %-8s\n', ...
           '(Hz)', 'Note', 'Note', 'Ratio', '(dBFS)', '(MIDI)', '');
    printf('%s\n', repmat('-', 1, 85));

    % Print each harmonic
    for i = 1:length(peak_freqs)
        freq = peak_freqs(i);
        amp = peak_amps(i);
        ratio = freq / fundamental;
        note_name = freq_to_midi_note(freq);

        % Calculate MIDI note number
        midi_note = round(69 + 12 * log2(freq / 440));

        velocity = amplitude_to_midi(amp);
        dynamic = midi_to_dynamic(velocity);

        printf('%-12.1f %-10s %-10d %-13.2f %-12.1f %-10d %-8s\n', ...
               freq, note_name, midi_note, ratio, amp, velocity, dynamic);
    end
end

function plot_spectrum(signal, fs, window_size, title_str, freq_min, freq_max)
    if nargin < 3
        window_size = 8192;
    end
    if nargin < 4
        title_str = 'Spectrum';
    end
    if nargin < 5
        freq_min = 100;
    end
    if nargin < 6
        freq_max = 6000;
    end

    % Take chunk from middle of the file
    mid_point = floor(length(signal) / 2);
    start_idx = mid_point - floor(window_size / 2);
    chunk = signal(start_idx:start_idx + window_size - 1) .* hanning(window_size);

    spectrum = abs(fft(chunk));
    spectrum = spectrum / (window_size / 2);  % Normalize for dBFS
    spectrum_db = 20 * log10(spectrum(1:window_size/2) + 1e-10);
    freq_axis = (0:window_size/2-1)' * fs / window_size;

    plot(freq_axis, spectrum_db);
    xlim([freq_min freq_max]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dBFS)');
    title(title_str);
    grid on;
end

function plot_spectrogram(signal, fs, title_str, freq_min, freq_max)
    if nargin < 3
        title_str = 'Spectrogram';
    end
    if nargin < 4
        freq_min = 100;
    end
    if nargin < 5
        freq_max = 6000;
    end

    window = 1024;
    overlap = 768;
    nfft = 2048;

    [S, f, t] = specgram(signal, nfft, fs, window, overlap);

    % Normalize for dBFS
    S_normalized = abs(S) / (nfft / 2);
    S_db = 20 * log10(S_normalized + 1e-10);

    imagesc(t, f, S_db);
    axis xy;
    ylim([freq_min freq_max]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(title_str);
    c = colorbar;
    ylabel(c, 'Magnitude (dBFS)');
end

printf('Static frequency analysis functions loaded!\n');
