function diagnostics = waveform_demo()
%WAVEFORM_DEMO Draw RRC-shaped Gray-QPSK samples, real RF, FFT and eyes.
% Run waveform_demo(); after setting all parameters in params.m.
% Uses no toolbox: manual Hann, fft, plotting primitives, RRC FIR defined
% analytically in rrc_filter.m.
%
% Fourier interpretation:
%   Convolution with RRC pulse in the time domain multiplies the spectrum
%   by the RRC frequency response; its square (Tx*Rx) is raised cosine.
%   Complex mixing with exp(j*2*pi*fc*t) translates baseband around fc.
%   FFT of a Hann-windowed finite capture displays this spectral shift and
%   approximate occupied RF bandwidth Rs*(1+alpha), NOT an ideal brickwall.

p=params();
rng(p.rng_seed,'twister');
numBits=2*p.demo_symbols;
txBits=randi([0,1],numBits,1);
txSymbols=qpsk_mod(txBits);
[rxSymbols,diagnostics]=qpsk_passband_link(txSymbols,p.demo_EbN0_dB,p,true);
rxBits=qpsk_demod(rxSymbols);
errors=sum(txBits~=rxBits);
Rs=p.fs_Hz/p.sps;
B=Rs*(1+p.rrc_alpha);
if ~exist(p.figure_dir,'dir')
    mkdir(p.figure_dir);
end

% (a) 20 symbol periods of pulse-shaped complex baseband I(t), Q(t).
N=min(p.plot_baseband_symbols,p.demo_symbols)*p.sps;
t_ms=(0:N-1)'/p.fs_Hz*1e3;
fig=figure('Color','w','Position',[100 100 920 500]);
plot(t_ms,real(diagnostics.tx_bb_samples(1:N)),'-','LineWidth',1.2);
hold on;
plot(t_ms,imag(diagnostics.tx_bb_samples(1:N)),'-','LineWidth',1.2);
grid on;
xlabel('Time (ms)');
ylabel('Baseband amplitude (normalized)');
title(sprintf('RRC-Shaped QPSK Complex Baseband: First %d Symbols',N/p.sps));
legend('In-phase I(t)','Quadrature Q(t)','Location','best');
pathA=fullfile(p.figure_dir,p.figure_bb_time);
print(fig,pathA,'-dpng',sprintf('-r%d',p.png_dpi));

% (b) 5 symbol periods of *real* 40 kHz RF waveform.
Nrf=min(p.plot_passband_symbols,p.demo_symbols)*p.sps;
t_ms=(0:Nrf-1)'/p.fs_Hz*1e3;
fig=figure('Color','w','Position',[100 100 920 460]);
plot(t_ms,diagnostics.tx_passband_samples(1:Nrf),'-','LineWidth',1.4);
grid on;
xlabel('Time (ms)');
ylabel('Real passband amplitude (normalized)');
title(sprintf('Real QPSK RF: First %d Symbols, f_c = %.0f kHz', ...
    Nrf/p.sps,p.fc_Hz/1e3));
pathB=fullfile(p.figure_dir,p.figure_rf_time);
print(fig,pathB,'-dpng',sprintf('-r%d',p.png_dpi));

% (c) Hann-windowed complex baseband and real passband FFT spectra.
[f,baseDB]=hann_fft_db(diagnostics.tx_bb_samples,p.fs_Hz,p.spectrum_nfft);
[~,passDB]=hann_fft_db(diagnostics.tx_passband_samples,p.fs_Hz,p.spectrum_nfft);
fig=figure('Color','w','Position',[100 100 1000 700]);
subplot(2,1,1);
plot(f/1e3,baseDB,'LineWidth',1.0); hold on;
grid on; ylim([-110 5]); xlim([-2*B/1e3,2*B/1e3]);
markY=ylim;
plot([-B/2 -B/2]/1e3,markY,'--','LineWidth',1.0);
plot([ B/2  B/2]/1e3,markY,'--','LineWidth',1.0);
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB relative to peak)');
title(sprintf('Baseband FFT (Hann): Nominal Full Width %.2f kHz',B/1e3));
legend('Baseband spectrum','RRC band edges','Location','southwest');
subplot(2,1,2);
plot(f/1e3,passDB,'LineWidth',1.0); hold on;
grid on; ylim([-110 5]); xlim([0 p.fs_Hz/(2*1e3)]);
markY=ylim;
plot([p.fc_Hz p.fc_Hz]/1e3,markY,'--','LineWidth',1.25);
plot([p.fc_Hz-B/2 p.fc_Hz-B/2]/1e3,markY,':','LineWidth',1.2);
plot([p.fc_Hz+B/2 p.fc_Hz+B/2]/1e3,markY,':','LineWidth',1.2);
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB relative to peak)');
title(sprintf('Real RF FFT (Hann): f_c=%.1f kHz, Nominal B=%.2f kHz', ...
    p.fc_Hz/1e3,B/1e3));
legend('RF spectrum','Carrier f_c','Nominal band edges','Location','southwest');
pathC=fullfile(p.figure_dir,p.figure_fft);
print(fig,pathC,'-dpng',sprintf('-r%d',p.png_dpi));

% (d) Eye diagram AFTER receive-matched RRC, at selected Eb/N0.
% The two filters total group delay 2*gd; the link function already
% aligns this waveform to symbol #1. Decision instants lie at half and
% 1.5 of these eye intervals, giving maximum eye opening with zero ISI.
numTraces=min(p.eye_traces,p.demo_symbols-3);
span=2*p.sps;
xEye=(0:span)'/p.sps;
fig=figure('Color','w','Position',[100 100 1000 480]);
subplot(1,2,1); hold on;
for k=1:numTraces
    first=1+p.sps/2+(k-1)*p.sps;
    plot(xEye,real(diagnostics.rx_matched_samples(first:first+span)), ...
        'LineWidth',0.45);
end
grid on;
xlabel('Time (symbol periods)');
ylabel('Matched-filter I amplitude (normalized)');
title('Receive Eye Diagram: In-phase');
subplot(1,2,2); hold on;
for k=1:numTraces
    first=1+p.sps/2+(k-1)*p.sps;
    plot(xEye,imag(diagnostics.rx_matched_samples(first:first+span)), ...
        'LineWidth',0.45);
end
grid on;
xlabel('Time (symbol periods)');
ylabel('Matched-filter Q amplitude (normalized)');
title(sprintf('Receive Eye Diagram: Quadrature, E_b/N_0=%g dB', ...
    p.demo_EbN0_dB));
pathD=fullfile(p.figure_dir,p.figure_eye);
print(fig,pathD,'-dpng',sprintf('-r%d',p.png_dpi));

% Numerical summary of the waveform experiment (not only figures).
fprintf('\nROBOT NURSE | RRC-RF DIAGNOSTIC SUMMARY\n');
fprintf('%-32s | %14s | %-12s\n','Quantity','Value','Unit');
fprintf('%s\n',repmat('-',1,66));
fprintf('%-32s | %14.2f | %-12s\n','Sampling frequency',p.fs_Hz,'Hz');
fprintf('%-32s | %14d | %-12s\n','Samples per symbol',p.sps,'samples/sym');
fprintf('%-32s | %14.2f | %-12s\n','Symbol rate',Rs,'sym/s');
fprintf('%-32s | %14.2f | %-12s\n','QPSK bit rate',Rs*p.bits_per_sym,'bit/s');
fprintf('%-32s | %14.2f | %-12s\n','Carrier',p.fc_Hz,'Hz');
fprintf('%-32s | %14.2f | %-12s\n','Nominal RF occupied width',B,'Hz');
fprintf('%-32s | %14.2f | %-12s\n','Tx RRC group delay',diagnostics.filter_group_delay_samples,'samples');
fprintf('%-32s | %14.2f | %-12s\n','Tx+Rx total group delay',diagnostics.total_delay_samples,'samples');
fprintf('%-32s | %14.2f | %-12s\n','Demo Eb/N0',p.demo_EbN0_dB,'dB');
fprintf('%-32s | %14.2f | %-12s\n','RF SNR per sample',diagnostics.sample_snr_dB,'dB');
fprintf('%-32s | %14d | %-12s\n','Demo bit errors',errors,'bits');
fprintf('%-32s | %14.6g | %-12s\n','Demo measured BER',errors/numBits,'dimensionless');
fprintf('Saved 4 figures at %d dpi in %s\n',p.png_dpi,p.figure_dir);

diagnostics.demo_bit_errors=errors;
diagnostics.demo_bits=numBits;
diagnostics.demo_BER=errors/numBits;
diagnostics.paths={pathA,pathB,pathC,pathD};
end
