function [frequency_Hz, amplitude_dB] = hann_fft_db(x, fs_Hz, nfft)
%HANN_FFT_DB Centered magnitude spectrum with manually implemented Hann.
%   [f,A] = hann_fft_db(x,fs,nfft)
%   A = 20log10(|FFT(x.*w)|/max(|FFT(x.*w)|)), in dB relative to peak.
%   The Hann window w[n]=0.5-0.5*cos(2*pi*n/(N-1)) reduces spectral leakage.
%   FFT (Discrete Fourier Transform) visualizes pulse-shaping bandwidth,
%   carrier shifting and the double-sided RF spectrum WITHOUT pwelch.
% No Signal Processing Toolbox functions (e.g. hann, hanning, pwelch).

if ~isvector(x) || numel(x)<2 || ~(isscalar(fs_Hz) && fs_Hz>0)
    error('hann_fft_db:Input','x must have >=2 samples and fs must be positive.');
end
if ~(isscalar(nfft) && nfft>=numel(x) && mod(nfft,2)==0)
    error('hann_fft_db:Nfft','nfft must be even and >= numel(x).');
end
x = x(:);
N = numel(x);
n=(0:N-1).';
w=0.5-0.5*cos(2*pi*n/(N-1));
X=fftshift(fft(x.*w,nfft));
magnitude=abs(X);
peak=max(magnitude);
amplitude_dB=20*log10(max(magnitude/max(peak,eps),1e-12));
frequency_Hz=(-nfft/2:nfft/2-1).'*fs_Hz/nfft;
end
