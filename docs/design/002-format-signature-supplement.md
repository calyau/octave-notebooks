# Supplemental Design: Format Signature for Spectral .mat Files

## Overview

Add a format signature variable to all saved .mat files to enable reliable detection by the `mat2sdif` CLI tool. This ensures the general-purpose tool can confidently identify spectral analysis files and apply appropriate conversion logic.

## Rationale

The `mat2sdif` tool is designed as a **general-purpose .mat to SDIF converter** that:
1. Works with any MATLAB v5 .mat file containing numeric arrays
2. Auto-detects spectral analysis files from this DSP package
3. Applies smart defaults for spectral → SDIF conversion

The format signature provides high-confidence detection while maintaining the general-purpose nature of the tool.

## Changes Required

### 1. Update `data/packages/custom-dsp/spectral-io.m`

#### In `save_static_spectrum()` function:

**Add before the `save()` call:**

```octave
    % Format signature for mat2sdif tool detection
    octave_dsp_format = 'spectral_static';
    octave_dsp_version = 1.0;
```

**Update the `save()` call to include these variables:**

```octave
    save('-v6', filename, 'freqs', 'amps', 'num_partials', ...
         'fundamental', 'sample_rate', 'time_stamp', ...
         'instrument_id', 'midi_note', ...
         'octave_dsp_format', 'octave_dsp_version');
```

#### In `save_spectral_evolution()` function:

**Add before the `save()` call:**

```octave
    % Format signature for mat2sdif tool detection
    octave_dsp_format = 'spectral_evolution';
    octave_dsp_version = 1.0;
```

**Update the `save()` call to include these variables:**

```octave
    save('-v6', filename, 'times', 'freqs', 'amps', ...
         'num_frames', 'num_partials_per_frame', ...
         'fundamental', 'sample_rate', 'window_size', 'hop_size', ...
         'instrument_id', 'midi_note', ...
         'octave_dsp_format', 'octave_dsp_version');
```

## Format Signature Specification

### Variables Added

**`octave_dsp_format`** (string stored as ASCII array)
- Type: Character array (will be stored as numeric array of ASCII values)
- Possible values:
  - `'spectral_static'` - Single-frame spectral snapshot
  - `'spectral_evolution'` - Time-varying spectral analysis
- Purpose: Identifies the specific data structure/schema

**`octave_dsp_version`** (numeric)
- Type: Scalar double
- Value: `1.0` (current version)
- Purpose: Allows future format evolution with backward compatibility

### How They're Stored

When Octave saves a string variable to a .mat file, it's stored as a numeric array of ASCII values:

```octave
% In Octave:
octave_dsp_format = 'spectral_static';

% In .mat file (as seen by matfile crate):
% Array of doubles: [115, 112, 101, 99, 116, 114, 97, 108, ...]
```

The Rust `matfile` crate will read this as `Array<f64>`, which can be converted back to a string:

```rust
let format_bytes = mat.find_by_name("octave_dsp_format")?;
let format_str = String::from_utf8(
    format_bytes.iter().map(|&x| x as u8).collect()
)?;
```

## Detection Logic for mat2sdif

The Rust tool will use this two-stage detection:

### Stage 1: Check for format signature (high confidence)

```rust
fn detect_octave_dsp_format(mat: &MatFile) -> Option<OctaveDspFormat> {
    // Try to read format signature
    if let Ok(format_array) = mat.find_by_name("octave_dsp_format") {
        // Convert numeric array to string
        let format_str = String::from_utf8(
            format_array.iter().map(|&x| x as u8).collect()
        ).ok()?;
        
        // Read version
        let version = mat.find_by_name("octave_dsp_version")
            .ok()?.get(0).copied()?;
        
        match format_str.as_str() {
            "spectral_static" => Some(OctaveDspFormat::SpectralStatic { version }),
            "spectral_evolution" => Some(OctaveDspFormat::SpectralEvolution { version }),
            _ => None
        }
    } else {
        None
    }
}
```

### Stage 2: Fallback heuristic detection (for files without signature)

If no signature is found, fall back to checking for the variable combination:

```rust
fn detect_by_variables(mat: &MatFile) -> Option<SpectralFormat> {
    let has_freqs = mat.find_by_name("freqs").is_ok();
    let has_amps = mat.find_by_name("amps").is_ok();
    let has_fundamental = mat.find_by_name("fundamental").is_ok();
    let has_sample_rate = mat.find_by_name("sample_rate").is_ok();
    
    if has_freqs && has_amps && has_fundamental && has_sample_rate {
        if mat.find_by_name("times").is_ok() {
            Some(SpectralFormat::Evolution)
        } else if mat.find_by_name("num_partials").is_ok() {
            Some(SpectralFormat::Static)
        } else {
            None
        }
    } else {
        None
    }
}
```

## Updated File Format Specification

### Static Spectrum .mat Variables

| Variable | Type | Dimensions | Description |
|----------|------|------------|-------------|
| `freqs` | double | 1×N | Frequencies in Hz |
| `amps` | double | 1×N | Amplitudes in dBFS |
| `num_partials` | double | 1×1 | Number of partials |
| `fundamental` | double | 1×1 | Fundamental frequency (Hz) |
| `sample_rate` | double | 1×1 | Sample rate (Hz) |
| `time_stamp` | double | 1×1 | Time in seconds |
| `instrument_id` | double | 1×1 | Instrument code (0-255) |
| `midi_note` | double | 1×1 | MIDI note number (0-127) |
| `octave_dsp_format` | double | 1×K | ASCII: 'spectral_static' |
| `octave_dsp_version` | double | 1×1 | Format version (1.0) |

### Evolution Spectrum .mat Variables

| Variable | Type | Dimensions | Description |
|----------|------|------------|-------------|
| `times` | double | 1×M | Time points (seconds) |
| `freqs` | double | M×K | Frequencies (Hz), zero-padded |
| `amps` | double | M×K | Amplitudes (dBFS), zero-padded |
| `num_frames` | double | 1×1 | Number of frames |
| `num_partials_per_frame` | double | 1×M | Partials per frame |
| `fundamental` | double | 1×1 | Fundamental frequency (Hz) |
| `sample_rate` | double | 1×1 | Sample rate (Hz) |
| `window_size` | double | 1×1 | FFT window size (samples) |
| `hop_size` | double | 1×1 | Hop size (samples) |
| `instrument_id` | double | 1×1 | Instrument code (0-255) |
| `midi_note` | double | 1×1 | MIDI note number (0-127) |
| `octave_dsp_format` | double | 1×K | ASCII: 'spectral_evolution' |
| `octave_dsp_version` | double | 1×1 | Format version (1.0) |

## Benefits

1. **Reliable detection**: No ambiguity about file format
2. **Versioning support**: Can evolve format while maintaining backward compatibility
3. **General-purpose tool**: `mat2sdif` can safely fall back to generic conversion for non-spectral files
4. **Low overhead**: Two small variables add minimal file size
5. **Future-proof**: Easy to add new format types (e.g., `'partial_tracking'`, `'multi_instrument'`)

## Testing

After implementation, verify with:

```octave
% Create and save test file
meta.fundamental = 440.0;
meta.sample_rate = 48000;
save_static_spectrum('test.mat', [440, 880], [30, 25], meta);

% Verify signature variables exist
whos -file test.mat

% Expected output should include:
%   octave_dsp_format    1xK     double
%   octave_dsp_version   1x1     double
```

## Implementation Checklist

- [ ] Update `save_static_spectrum()` to include signature variables
- [ ] Update `save_spectral_evolution()` to include signature variables
- [ ] Test static spectrum save/load
- [ ] Test evolution spectrum save/load
- [ ] Verify signature variables are readable as ASCII arrays
- [ ] Update example scripts if needed
- [ ] Document in main README

## Future Format Types

If additional spectrum types are added later, use this naming convention:

- `'spectral_static'` - Single snapshot (current)
- `'spectral_evolution'` - Time-varying (current)
- `'partial_tracking'` - Individual sinusoid tracking
- `'multi_instrument'` - Multiple instruments analyzed together
- `'spectral_morph'` - Morphing between states

Version number should increment for breaking changes to existing formats:
- `1.0` - Initial version
- `1.1` - Backward-compatible additions
- `2.0` - Breaking changes

---

## Notes for mat2sdif Implementation

When implementing the Rust side:

1. **Always check signature first** - it's the most reliable method
2. **Fall back to heuristics** - for backward compatibility with files that don't have signatures
3. **Validate schema** - even with signature, verify expected variables exist
4. **Version checking** - warn if version is newer than supported
5. **Graceful degradation** - offer generic conversion mode if spectral detection fails

Example error messages:

```
✓ Detected Octave DSP spectral_evolution format (v1.0)
✓ Converting 145 frames to SDIF...

⚠ Unknown Octave DSP version 2.5 (supported: 1.x)
  Attempting conversion anyway...

✗ Spectral format detected but missing required variable: 'fundamental'
  Try: mat2sdif --generic input.mat output.sdif --help

ℹ No spectral format detected
  Use --help to see generic conversion options
```
