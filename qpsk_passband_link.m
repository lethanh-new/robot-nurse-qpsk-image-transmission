function [rxSymbols, diag] = qpsk_passband_link(txSymbols, EbN0_dB, p, add_noise)
%QPSK_PASSBAND_LINK Full sampled Gray-QPSK link, coherent real RF + AWGN.
%   [rxSymbols, diag] = qpsk_passband_link(txSymbols, EbN0_dB, p, true)
% Inputs are one complex QPSK symbol per bit pair (Es=1). Output contains
% one recovered matched-filter symbol per input symbol; call qpsk_demod.
%
% Tx: symbol -> zero-insert upsample (sps=16) -> unit-energy RRC FIR ->
%     real RF s[n] = sqrt(2)*Re{xBB[n]*exp(j*2*pi*fc*n/fs)}.
% Ch: sRF + real AWGN.  Rx: sqrt(2)*rRF*exp(-j*2*pi*fc*n/fs)
%     -> matched RRC FIR -> sample after TOTAL Tx+Rx group delay.
%
% ENERGY / NOISE SCALING (not MATLAB awgn()):
%   Es=1, Eb=0.5, gamma=10^(EbN0_dB/10), N0=Eb/gamma.
%   For unit-energy RRC, real RF noise variance PER SAMPLE is N0/2.
%   sqrt(2) downconversion followed by unit-energy matched RRC makes the
%   COMPLEX noise variance at the decision time approximately N0.
%   RF sample SNR = 2*Es/(sps*N0) = 4*gamma/sps.
%   In dB: SNR_RF_sample = EbN0_dB + 10*log10(4/sps).
%   (Complex-BB sample SNR = EbN0_dB + 10*log10(2/sps).)
% Carrier sync, symbol timing, and channel gain are assumed perfectly known.
%
% FFT is used here only as a fast, EXACT linear-FIR implementation:
% zero-padding to >= length(signal)+length(h)-1 avoids circular aliasing.
% The Fourier transform converts convolution in time to multiplication in
% frequency; waveform_demo.m also uses FFT to measure spectral occupancy.

if nargin < 4
    add_noise = true;
end
if ~isnumeric(txSymbols) || ~isvector(txSymbols) || isempty(txSymbols)
    error('qpsk_passband_link:Input', 'txSymbols must be a nonempty vector.');
end
if ~(isscalar(EbN0_dB) && isfinite(EbN0_dB))
    error('qpsk_passband_link:EbN0', 'EbN0_dB must be finite scalar.');
end
if p.Es ~= 1 || p.Eb ~= 0.5 || p.bits_per_sym ~= 2
    error('qpsk_passband_link:Energy', 'Expected Es=1, Eb=0.5, QPSK.');
end
if p.guard_symbols < p.rrc_span_symbols
    error('qpsk_passband_link:Guards', ...
        'guard_symbols must be >= rrc_span_symbols to discard transients.');
end
if p.channel_gain ~= 1
    error('qpsk_passband_link:Gain', ...
        'BER theory requires channel_gain=1 (AWGN-only baseline).');
end
if p.fc_Hz <= 0 || p.fc_Hz >= p.fs_Hz/2 || p.fs_Hz/p.sps <= 0
    error('qpsk_passband_link:Sampling', 'Invalid carrier/sampling frequencies.');
end
B = (p.fs_Hz/p.sps)*(1+p.rrc_alpha);
if p.fc_Hz-B/2 <= 0 || p.fc_Hz+B/2 >= p.fs_Hz/2
    error('qpsk_passband_link:Band', ...
        'RRC carrier band must fit inside the sampled positive RF band.');
end

symbols = txSymbols(:);
numData = numel(symbols);
[h, gd] = rrc_filter(p.rrc_alpha, p.rrc_span_symbols, p.sps);
G = p.guard_symbols;
% Real QPSK guard symbols avoid startup/ending transients in useful data.
guardSym = (1+1j)/sqrt(2);
framedSymbols = [repmat(guardSym,G,1); symbols; repmat(guardSym,G,1)];
numSamples = numel(framedSymbols)*p.sps;
up = zeros(numSamples,1);
up(1:p.sps:end) = framedSymbols;

% Two FFT-based, zero-padded FIR convolutions: no Signal Processing Toolbox.
nfft = 2^nextpow2(numSamples + numel(h) - 1);
H = fft(h,nfft);
filteredTx = ifft(fft(up,nfft).*H);
txBB = filteredTx(1:numSamples);

% Carrier starts at n=0 each frame. Each frame is an integer number of
% symbol periods; fc/Rs = 4 guarantees the same phase at every symbol.
n = (0:numSamples-1).';
if p.fc_Hz*4 == p.fs_Hz
    % At fc=fs/4, the carrier has the EXACT repeating phases 1,j,-1,-j;
    % this avoids trig roundoff and speeds up very long BER simulations.
    phases = [1;1j;-1;-1j];
    osc = phases(mod(n,4)+1);
else
    osc = exp(1j*2*pi*(p.fc_Hz/p.fs_Hz)*n);
end
sRF = sqrt(2)*real(txBB.*osc);

gamma_b = 10^(EbN0_dB/10);
N0 = p.Eb/gamma_b;
if add_noise
    rfNoise = sqrt(N0/2)*randn(numSamples,1); % REAL passband AWGN.
else
    rfNoise = zeros(numSamples,1);
end
rRF = sRF + rfNoise;
rxMixed = sqrt(2)*rRF.*conj(osc); % I/Q coherent downconverter.

filteredRx = ifft(fft(rxMixed,nfft).*H);
rxMF = filteredRx(1:numSamples);

% Signal path contains two causal RRC delays (one at Tx and one at Rx).
% MATLAB uses 1-based samples; first input symbol is placed at index 1.
dataSymStart = G*p.sps+1;
sampleIndices = dataSymStart + 2*gd + (0:numData-1)'*p.sps;
if sampleIndices(end) > numSamples
    error('qpsk_passband_link:Delay', 'Insufficient guard samples.');
end
rxSymbols = rxMF(sampleIndices);

if nargout > 1
    % Align diagnostic waveforms to data symbol #1; keep RF carrier and
    % TX pulse shapes at the physical TX group delay gd, not at 2*gd.
    samplesBB = dataSymStart + gd + (0:numData*p.sps-1).';
    samplesRX = dataSymStart + 2*gd + (0:numData*p.sps-1).';
    diag = struct();
    diag.tx_bb_samples = txBB(samplesBB);
    diag.tx_passband_samples = sRF(samplesBB);
    diag.rx_matched_samples = rxMF(samplesRX);
    diag.rx_rf_samples = rRF(samplesBB);
    diag.rrc_coeffs = h;
    diag.filter_group_delay_samples = gd;
    diag.total_delay_samples = 2*gd;
    diag.symbol_rate_Hz = p.fs_Hz/p.sps;
    diag.nominal_rf_bandwidth_Hz = (p.fs_Hz/p.sps)*(1+p.rrc_alpha);
    diag.N0 = N0;
    diag.rf_noise_variance = N0/2;
    diag.sample_snr_dB = EbN0_dB + 10*log10(4/p.sps);
else
    diag = [];
end
end
