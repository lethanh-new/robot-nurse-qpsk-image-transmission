function test = qpsk_chain_check()
%QPSK_CHAIN_CHECK Numerical sanity checks of the 16-sample/symbol RF link.
% Tests constellation labeling, RRC unit energy, timing delay, noiseless
% BER, and N0 recovery through real AWGN/coherent matched filter.
% This test is separate from the 100-error Monte Carlo BER experiment.
% In particular, no noise samples are tuned to force BER agreement.

p=params();
rng(p.rng_seed,'twister');
[h,gd]=rrc_filter(p.rrc_alpha,p.rrc_span_symbols,p.sps);
labels=[0;0; 0;1; 1;0; 1;1];
expected=[1+1j;1-1j;-1+1j;-1-1j]/sqrt(2);
modCheck=max(abs(qpsk_mod(labels)-expected));
demodCheck=sum(qpsk_demod(expected)~=labels);
errRRC=abs(sum(abs(h).^2)-1);

bits=randi([0 1],2*p.selftest_symbols,1);
symbols=qpsk_mod(bits);
[cleanSymbols,diag]=qpsk_passband_link(symbols,p.demo_EbN0_dB,p,false);
noiselessErrors=sum(qpsk_demod(cleanSymbols)~=bits);
[noisySymbols,~]=qpsk_passband_link(symbols,p.demo_EbN0_dB,p,true);
noiseDecision=noisySymbols-cleanSymbols;
noisePower=mean(abs(noiseDecision).^2);
noiseRelError=abs(noisePower/diag.N0-1);

fprintf('\nQPSK 16-SAMPLE RF LINK SELF CHECK\n');
fprintf('%-34s | %16s\n','Check','Measured value');
fprintf('%s\n',repmat('-',1,57));
fprintf('%-34s | %16.6g\n','Max Gray constellation error',modCheck);
fprintf('%-34s | %16d\n','Gray demap errors (4 symbols)',demodCheck);
fprintf('%-34s | %16.6g\n','RRC energy error from unity',errRRC);
fprintf('%-34s | %16d\n','TX group delay (samples)',gd);
fprintf('%-34s | %16d\n','TX+RX group delay (samples)',diag.total_delay_samples);
fprintf('%-34s | %16d\n','Noiseless payload bit errors',noiselessErrors);
fprintf('%-34s | %16.6g\n','Theoretical decision noise N0',diag.N0);
fprintf('%-34s | %16.6g\n','Measured decision noise E|n|^2',noisePower);
fprintf('%-34s | %15.3f%%\n','Noise power relative error',100*noiseRelError);

passed=(modCheck<1e-12) && (demodCheck==0) && (errRRC<1e-12) && ...
       (noiselessErrors==0) && (noiseRelError<=p.selftest_noise_tolerance);
if passed
    fprintf('PASS: consistency checks satisfied.\n');
else
    fprintf('CHECK FAILED: inspect the measured values above; do not alter data.\n');
end
test=struct('passed',passed,'noiseless_bit_errors',noiselessErrors, ...
    'rrc_energy_error',errRRC,'noise_power',noisePower, ...
    'expected_N0',diag.N0,'noise_rel_error',noiseRelError, ...
    'group_delay_samples',gd);
end
