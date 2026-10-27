# Robot Nurse -> Monitoring Station | QPSK with RRC + real RF

MATLAB only, no Signal Processing / Communications Toolbox or Simulink.

## Files

- `params.m` - all physical/model/Monte Carlo/plot parameters
- `qpsk_mod.m` - Gray mapping, Es=1 (1 complex symbol per 2 bits)
- `qpsk_demod.m` - hard decision with coherent timing and carrier
- `rrc_filter.m` - analytic root-raised cosine, including singularities; normalized energy
- `qpsk_passband_link.m` - oversample x16, TX RRC, real 40-kHz carrier,
  real passband AWGN, coherent downconvert, matched RRC, total group delay removal
- `hann_fft_db.m` - manual Hann + FFT spectrum (no `pwelch`, `hann`, `rcosdesign`)
- `waveform_demo.m` - 4 plots: 20 symbols I/Q, 5 passband symbols, spectra, eye
- `qpsk_chain_check.m` - quick numeric check of Gray mapping, unit RRC energy, timing and noise variance
- `ber_sim.m` - Monte Carlo full chain 0:1:10 dB, stop at >=100 errors or 2e7 bits

## Run

From the project folder in MATLAB:

```matlab
p = params();
test = qpsk_chain_check();
diagnostics = waveform_demo();
results = ber_sim();
```

The passband Monte Carlo at 10 dB can process **20 million bits**, or **160
million samples** before the bit cap; allow sufficient CPU time and RAM.
This run is intentionally more expensive than the prior 1-sample/symbol baseline.

To compare against the earlier baseline, edit only `params.m`:

```matlab
p.ber_mode = 'baseband_awgn';
```

Then run `ber_sim()` again. Restore `'passband_rrc'` for the complete chain.

## Signal and sample SNR conventions

Symbol energy Es=1, Eb=0.5, gamma=10^(EbN0/10), N0=Eb/gamma.

- TX RRC is **unit discrete-time energy**: sum(h.^2)=1.
- TX xBB is upsampled by inserting 15 zeros between symbols and RRC filtered.
- TX real passband: `sRF=sqrt(2)*real(xBB.*exp(1j*2*pi*fc*n/fs))`.
- RF AWGN: `sqrt(N0/2)*randn(...)` with variance **N0/2 per real sample**.
- RF sample power ~Es/sps, RF noise variance N0/2, hence
  `SNR_RF_sample_dB=EbN0_dB + 10*log10(4/sps)` = EbN0 - 6.0206 dB at sps=16.
- The downconverter multiplies by `sqrt(2)*exp(-1j*2*pi*fc*n/fs)` and
  RX matched filter produces approximately complex noise variance N0 at
  symbol decisions (the phase-conjugated 2fc image is suppressed by RRC).
- With perfectly aligned carrier, gain and timing, theoretical
  `Pb=0.5*erfc(sqrt(gamma))` should match Monte Carlo within random error.
- RRC span=10 symbols and total group delay = 160 samples.
- Guards prevent contamination by zero-padded filter startup and endpoint.

**No parameters or BER observations are changed to force agreement.** With a
100-error stop rule, a relative deviation beyond 20% is statistically
possible (~10% relative standard error), even for a physically correct link.
`ber_sim` prints whether every >=100-error point passes the tolerance audit;
points hitting the 20-million-bit cap below 100 errors are excluded from it.

## FFT note

FFT is the finite-sample discrete Fourier transform, used (1) to implement
exact linear convolutions with zero-padding as multiplication in frequency,
and (2) to display the translation of an RRC-limited complex baseband
spectrum to a real passband centered at +/-40 kHz. A manually generated
Hann window reduces spectral leakage for these plots. The nominal **full**
RF occupied width at alpha=0.35 is `Rs*(1+alpha)=13.5 kHz` around fc,
with edges 33.25 and 46.75 kHz on the positive-frequency side. Finite RRC
truncation and finite FFT capture cause sidelobes and do not change theory.

## Outputs

`figures/rrc_baseband_iq.png`, `figures/rrc_passband_5symbols.png`,
`figures/rrc_hann_fft_baseband_passband.png`, `figures/rrc_matched_eye.png`,
`figures/qpsk_ber_vs_theory.png`, `figures/qpsk_ber_table.csv`.
Figures use 200 dpi, labelled axes with units, descriptive titles.

## Image transmission stage (128x128 grayscale)

Prepared sample images are included under `images/`.
Default experiment image: `images/01_monitor_active_128.png`.

Run:

```matlab
clear; clc; close all;
results_img = image_tx_demo();
```

Optional alternate image:

```matlab
results_img = image_tx_demo(fullfile('images','06_ecg_flatline_128.png'));
```

This stage saves the following outputs in `figures/`:

- `image_reconstruction_grid.png` : one figure with original image and reconstructed images.
- `image_compare_*dB.png` : one comparison figure per Eb/N0 point.
- `image_restored_*dB.png` : raw reconstructed 128x128 PNG.
- `image_results_table.csv` : numerical table with BER, MSE, and PSNR.

The image BER comes from one finite image frame (131072 bits), so it can differ from the long-run Monte Carlo BER curve. Do not "adjust" results.

### MATLAB `bitget` compatibility correction

`image_to_bits.m` now extracts 8 bit planes using scalar bit positions
(`for k = 1:8`), avoiding the size mismatch raised by
`bitget(pixels, 8:-1:1)` on some MATLAB releases.

Run the lossless input/output sanity check before the RF experiment:

```matlab
assert(image_self_check());
results_img = image_tx_demo();
```

### PNG output directory fix (MATLAB on Windows)

`image_tx_demo` now resolves `figures/` from the location of `image_tx_demo.m`
using `fileparts(mfilename('fullpath'))`. It creates that absolute directory
and tests a small PNG write BEFORE running the RF simulation. If the test
fails, check Windows folder permissions or update `p.figure_dir` in
`params.m` to a writable directory name. This does not change any random,
modulation, SNR, or image settings.

## Synthetic 64x64 PCM8 vs JPEG encoding experiment (new `image_tx.m`)

This stage uses a **new entirely synthetic 64x64 8-bit grayscale graphic**:
there is no actual patient data, and it does not use the 128x128 example images.
It uses the **same full stage-2 QPSK/RRC/40-kHz real-RF + AWGN receiver**.
The fixed experiment levels are Eb/N0 = **0, 4, 8 dB**, with `rng(1)`.

```matlab
clear; clc; close all;
results = image_tx();
```

All physical and encoding settings are in `params.m`, including JPEG quality
(default = 75), `p.image64_ebn0_dB`, the image size and output file names.

Two payload types are sent independently:

1. **PCM8:** serialize the 64x64 `uint8` pixels, 8 bits MSB-first per pixel.
   Uncompressed frame size = **32768 bits**. Pixels are reconstructed from
   received bytes regardless of bit errors.
2. **JPEG:** `imwrite(...,'jpg','Quality',...)` creates an actual lossy JPEG.
   The **entire JPEG file** is read as bytes, transmitted bit-for-bit MSB
   first (including its header), then written back as bytes and decoded by
   `imread`. If the corrupted JPEG cannot be decoded, `try/catch` reports
   `DECODE_FAIL` and MSE/PSNR are `NaN`. If the decoder returns a warning,
   status is `DECODE_WARN`. A successful decode with bit errors is not
   guaranteed to look correct.

Payload time = `number_of_bits / (fs/sps * bits_per_sym)`. Here Rb = **20000
bit/s**, so the PCM8 payload time is **1.6384 s**. The JPEG payload size
and corresponding time are measured from the **actual MATLAB-created JPEG**;
they are **not hard-coded**. The duration excludes filter guard symbols,
packet headers and channel framing. The JPEG compression-only MSE and PSNR
are also reported because JPEG is lossy even with zero channel bit errors.

Outputs (in absolute project-local `figures/` directory):

- `synthetic_monitor_64x64.png` (raw generated input).
- `synthetic_monitor_q75.jpg` (source JPEG file).
- `image_tx_64_pcm_comparison.png`: original + reconstructions 0/4/8 dB.
- `image_tx_64_jpeg_comparison.png`: original + reconstructions or explicit
  `JPEG DECODE FAIL` panels 0/4/8 dB.
- `image_tx_64_results.csv`: 6 rows, input bits, errors, measured BER,
  MSE, PSNR, payload transmission time, number of bytes and decode status.
- Per-case recovered PNGs where decoding succeeds and the received JPEG
  byte files used by the decoder.

No EVM, AWGN, or BER values are altered to make the plots more attractive.
The one-frame BER can fluctuate relative to a long Monte Carlo benchmark.

New helper functions (one function per `.m`):
`make_synthetic_image64.m`, `bytes_to_bits.m`, `bits_to_bytes.m`.
