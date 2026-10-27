function p = params()
%PARAMS All tunable parameters for Robot Nurse QPSK experiments.
% MATLAB base functions only: no Communications/Signal Processing Toolbox.
% Unit-energy Gray QPSK, complex baseband and real passband modes.
%
% Frequency convention: 16 samples/symbol at fs=160 kHz gives Rs=10 ksym/s.
% With RRC roll-off 0.35, nominal RF occupied width = Rs*(1+alpha).

% Repeatable Monte Carlo, channel and modulation.
p.rng_seed        = 1;
p.Es              = 1;
p.bits_per_sym    = 2;
p.Eb              = p.Es / p.bits_per_sym;
p.EbN0_dB         = 0:1:10;
p.target_errors   = 100;
p.max_bits        = 2e7;
p.chunk_bits      = 100000;     % Even, bounded memory; dynamically shortened.
p.chunk_margin    = 1.5;        % Only optimizes frame length (never changes stop rule).
p.ber_mode        = 'passband_rrc'; % 'passband_rrc' or 'baseband_awgn'.
p.ber_tolerance_percent = 20;   % Report PASS/FAIL, never modify outcomes.

% Matched RRC pulse shaping, sampling and coherent RF link.
p.sps             = 16;         % Samples/symbol.
p.fs_Hz           = 160000;     % Sample rate [Hz].
p.fc_Hz           = 40000;      % Real RF carrier [Hz].
p.rrc_alpha       = 0.35;       % Root-raised-cosine roll-off, 0 <= a <= 1.
p.rrc_span_symbols = 10;        % Even; total RRC filter length = span*sps + 1.
p.guard_symbols   = 10;         % >= RRC span so frames are independent at decisions.
p.channel_gain    = 1;         % Ideal AWGN benchmark: known carrier, no fading.

% Plot/demo parameters; demo bits are distinct from BER experiment.
p.demo_symbols    = 1200;
p.selftest_symbols = 4000;     % QPSK symbols used for finite-noise-power test.
p.selftest_noise_tolerance = 0.15; % Check measured decision noise energy.
p.demo_EbN0_dB    = 8;
p.plot_baseband_symbols = 20;
p.plot_passband_symbols = 5;
p.eye_traces      = 160;
p.spectrum_nfft   = 65536;
p.figure_dir      = 'figures';
p.figure_name     = 'qpsk_ber_vs_theory.png';
p.figure_bb_time  = 'rrc_baseband_iq.png';
p.figure_rf_time  = 'rrc_passband_5symbols.png';
p.figure_fft      = 'rrc_hann_fft_baseband_passband.png';
p.figure_eye      = 'rrc_matched_eye.png';
p.png_dpi         = 200;
% Synthetic 64x64 image experiment (no patient images or private data).
p.image64_side       = 64;                 % Pixel width and height.
p.image64_ebn0_dB    = [0 4 8];           % Exactly three Eb/N0 test points.
p.jpeg_quality       = 75;                % JPEG lossy quality [0,100].
p.image64_png        = 'synthetic_monitor_64x64.png';
p.image64_jpeg       = 'synthetic_monitor_q75.jpg';
p.image64_csv        = 'image_tx_64_results.csv';
p.image64_pcm_fig    = 'image_tx_64_pcm_comparison.png';
p.image64_jpeg_fig   = 'image_tx_64_jpeg_comparison.png';
p.image64_recovered_prefix = 'image_tx_64_recovered_';
p.image64_rx_jpeg_prefix = 'image_tx_64_rx_';
end
