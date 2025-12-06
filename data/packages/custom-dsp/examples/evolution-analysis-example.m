% Example: Time-varying spectrum analysis
pkg load signal;

% Load custom DSP functions
addpath('../');
run('../freq-analysis.m');
run('../spectral-io.m');
run('../time-varying-analysis.m');

% Analyze oboe A4 evolution
printf('Analyzing oboe evolution...\n');
[times, freqs, amps, meta] = analyze_spectral_evolution('../../../audio/oboe.wav', 4096, 2048, 15);

% Add instrument metadata
meta.fundamental = 439.5;
meta.instrument_id = 1;  % oboe
meta.midi_note = 69;     % A4

% Save to .mat file
save_spectral_evolution('../../../output/oboe_A4_evolution.mat', times, freqs, amps, meta);

printf('\nAnalysis complete!\n');
printf('Frames: %d\n', meta.num_frames);
printf('Duration: %.2f seconds\n', times(end) - times(1));
printf('Time resolution: %.3f seconds\n', meta.hop_size / meta.sample_rate);
