% Example: Static spectrum analysis
pkg load signal;

% Load custom DSP functions
addpath('../');
run('../freq-analysis.m');
run('../spectral-io.m');
run('../time-varying-analysis.m');

% Analyze oboe A4
[freqs, amps] = analyze_static_spectrum('../../../audio/oboe.wav', 8192);

% Display results
printf('\nOboe A4 Static Analysis:\n');
print_harmonics(freqs, amps);

% Prepare metadata
meta.fundamental = 439.5;
meta.sample_rate = 48000;
meta.instrument_id = 1;  % 1 = oboe
meta.midi_note = 69;     % A4

% Save to .mat file
save_static_spectrum('../../../output/oboe_A4_static.mat', freqs, amps, meta);
