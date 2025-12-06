# Octave Jupyter DSP Environment

A Docker-based Jupyter notebook environment with Octave kernel and signal processing packages.

## Quick Start

```bash
docker compose up --build
```

Then open <http://localhost:8889> in your browser.

**Note:** This runs on port 8889 so you can run it alongside the Rust Jupyter setup on 8888.

## Included Packages

- **signal** - DSP functions (fft, filters, windows, spectrogram)
- **image** - Image processing
- **statistics** - Statistical functions
- **control** - Control systems toolbox
- **io** - File I/O utilities

## Key Functions for Harmonic Analysis

| Function | Description |
|----------|-------------|
| `audioread(file)` | Load WAV/FLAC/OGG files |
| `fft(x)` | Fast Fourier Transform |
| `specgram(x)` | Compute spectrogram |
| `findpeaks(x)` | Find local maxima |
| `hanning(n)` | Hanning window |
| `butter(n, wc)` | Butterworth filter design |

## Loading Your Audio Files

Put your WAV files in the `data/` folder, then:

```octave
[y, fs] = audioread('clarinet.wav');

% If stereo, take first channel
if size(y, 2) > 1
    y = y(:, 1);
end

% Analyze harmonics
[freqs, amps] = find_harmonics(y, fs, 8192);
print_harmonics(freqs, amps);
```

## Stopping

```bash
docker compose down
```
