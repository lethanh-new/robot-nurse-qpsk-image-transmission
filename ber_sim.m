function results = ber_sim()
%BER_SIM QPSK Monte Carlo BER through FULL RRC -> real RF -> AWGN -> RRC.
% Run: results = ber_sim();  (or change p.ber_mode to 'baseband_awgn').
% Does not need Communications Toolbox or Signal Processing Toolbox.
%
% Symbol Es=1, Eb=0.5, N0=Eb/10^(EbN0_dB/10). Gray QPSK.
% For the full real-passband sampled system, unit-energy RRC gives:
%     var(n_RF[n])=N0/2, n_RF REAL; sampled SNR_RF=4*gamma_b/sps.
%     SNR_RF(dB)=EbN0(dB)+10*log10(4/sps), sps=16.
%   Coherent sqrt(2) downmix + unit-energy matched RRC returns output
%   COMPLEX noise approximately CN(0,N0), i.e. the same as baseband AWGN.
%   Compare Pb=Q(sqrt(2*gamma_b))=0.5*erfc(sqrt(gamma_b)).
%
% Stop exactly when the target-th bit error is reached in a frame OR once
% max_bits is reached. The Monte Carlo error count is NOT adjusted to make
% the graph match theory. A +/-20% check is only reported, not enforced.

p = params();
assert(p.Es == 1 && p.bits_per_sym == 2 && p.Eb == 0.5, ...
    'Expected Es=1 and Eb=0.5.');
assert(mod(p.chunk_bits,2)==0 && mod(p.max_bits,2)==0, ...
    'chunk_bits and max_bits must be even.');
assert(p.chunk_bits > 0 && p.max_bits > 0 && p.target_errors >= 100, ...
    'Invalid Monte Carlo limits.');
assert(any(strcmp(p.ber_mode,{'baseband_awgn','passband_rrc'})), ...
    'Unknown p.ber_mode.');
rng(p.rng_seed,'twister');

snr_dB = p.EbN0_dB(:);
points = numel(snr_dB);
numBits = zeros(points,1);
numErrors = zeros(points,1);
BER_sim = zeros(points,1);
BER_theory = zeros(points,1);
delta_pct = zeros(points,1);

fprintf('\nROBOT NURSE -> MONITORING STATION | QPSK BER\n');
fprintf('Mode: %s | rng(%d) | Es=%.1f | Eb=%.1f\n', ...
    p.ber_mode, p.rng_seed, p.Es, p.Eb);
if strcmp(p.ber_mode, 'passband_rrc')
    fprintf('fs=%.0f Hz | fc=%.0f Hz | sps=%d | alpha=%.3f | RRC span=%d symbols\n', ...
        p.fs_Hz,p.fc_Hz,p.sps,p.rrc_alpha,p.rrc_span_symbols);
    fprintf('Sample SNR_RF = Eb/N0(dB) + 10log10(4/sps) = Eb/N0(dB) %+.3f dB\n', ...
        10*log10(4/p.sps));
end
fprintf('Stop at >=%d errors or >=%.0f bits (per Eb/N0 point).\n\n', ...
    p.target_errors, p.max_bits);

for k = 1:points
    gamma = 10^(snr_dB(k)/10);
    N0 = p.Eb/gamma;
    BER_theory(k) = 0.5*erfc(sqrt(gamma));  % Q(x)=erfc(x/sqrt(2))/2.

    while numErrors(k) < p.target_errors && numBits(k) < p.max_bits
        remaining = p.target_errors - numErrors(k);
        expectedBits = 2*ceil(p.chunk_margin*remaining/(2*BER_theory(k)));
        blockLen = min([p.chunk_bits, p.max_bits-numBits(k), max(2,expectedBits)]);
        if mod(blockLen,2) ~= 0
            blockLen = blockLen - 1;
        end
        txBits = randi([0,1],blockLen,1);
        txSym = qpsk_mod(txBits);
        if strcmp(p.ber_mode,'passband_rrc')
            % The function performs ALL oversampled Tx/Rx operations.
            rxSym = qpsk_passband_link(txSym,snr_dB(k),p,true);
        else
            n = sqrt(N0/2)*(randn(size(txSym))+1j*randn(size(txSym)));
            rxSym = txSym+n;
        end
        rxBits = qpsk_demod(rxSym);
        errorLocations = find(txBits ~= rxBits);
        if numel(errorLocations) >= remaining
            % Include only received bits through the target-th error.
            numBits(k) = numBits(k) + errorLocations(remaining);
            numErrors(k) = p.target_errors;
        else
            numBits(k) = numBits(k) + blockLen;
            numErrors(k) = numErrors(k) + numel(errorLocations);
        end
    end
    BER_sim(k) = numErrors(k)/numBits(k);
    delta_pct(k) = 100*(BER_sim(k)-BER_theory(k))/BER_theory(k);
    fprintf('Processed Eb/N0=%2d dB: bits=%d, errors=%d, deviation=%+.2f%%\n', ...
        snr_dB(k),numBits(k),numErrors(k),delta_pct(k));
end

fprintf('\n%8s | %12s | %9s | %12s | %12s | %13s\n', ...
    'EbN0(dB)','Bits','Errors','BER_sim','BER_theory','Deviation(%)');
fprintf('%s\n', repmat('-',1,84));
for k=1:points
    fprintf('%8.1f | %12d | %9d | %12.6e | %12.6e | %+12.3f%%\n', ...
        snr_dB(k),numBits(k),numErrors(k),BER_sim(k),BER_theory(k),delta_pct(k));
end
pass_mask = numErrors >= p.target_errors;
pass_pct = abs(delta_pct) <= p.ber_tolerance_percent;
if any(pass_mask)
    fprintf('\nTolerance audit at >=%d errors: %d/%d within +/-%.1f%%.\n', ...
        p.target_errors, sum(pass_pct & pass_mask),sum(pass_mask), ...
        p.ber_tolerance_percent);
    if any(pass_mask & ~pass_pct)
        fprintf('WARNING: Some valid Monte Carlo points exceed tolerance.');
        fprintf(' Preserve measured results: 100 errors have statistical uncertainty.\n');
    end
end
if any(~pass_mask)
    fprintf('Cap reached with fewer than %d errors at Eb/N0:',p.target_errors);
    fprintf(' %g',snr_dB(~pass_mask));
    fprintf(' dB. These points are excluded from the >=100-error audit.\n');
end

if ~exist(p.figure_dir,'dir')
    mkdir(p.figure_dir);
end
fig = figure('Color','w','Position',[120 120 800 540]);
semilogy(snr_dB,BER_sim,'o','MarkerSize',7,'LineWidth',1.5);
hold on;
semilogy(snr_dB,BER_theory,'-','LineWidth',1.6);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit error rate, BER (dimensionless)');
if strcmp(p.ber_mode,'passband_rrc')
    title('Gray QPSK: RRC + 40 kHz Passband + AWGN + Matched RRC');
else
    title('Gray QPSK: 1 Sample/Symbol + Complex AWGN');
end
legend('Monte Carlo (markers)','Theory Q(sqrt(2E_b/N_0))', ...
    'Location','southwest');
outputPath = fullfile(p.figure_dir,p.figure_name);
print(fig,outputPath,'-dpng',sprintf('-r%d',p.png_dpi));

% Export raw numeric evidence for the Assignment report.
csvPath = fullfile(p.figure_dir,'qpsk_ber_table.csv');
fid = fopen(csvPath,'w');
if fid < 0
    error('ber_sim:CSV', 'Cannot open output CSV file.');
end
fprintf(fid,'EbN0_dB,num_bits,num_errors,BER_sim,BER_theory,deviation_percent\n');
for k=1:points
    fprintf(fid,'%.0f,%.0f,%.0f,%.12g,%.12g,%+.8f\n', ...
        snr_dB(k),numBits(k),numErrors(k),BER_sim(k),BER_theory(k),delta_pct(k));
end
fclose(fid);
fprintf('Saved: %s (%d dpi)\n',outputPath,p.png_dpi);
fprintf('Saved: %s\n',csvPath);

results = struct('EbN0_dB',snr_dB,'numBits',numBits, ...
    'numErrors',numErrors,'BER_sim',BER_sim,'BER_theory',BER_theory, ...
    'deviation_percent',delta_pct,'within_20pct',pass_pct, ...
    'has_100_errors',pass_mask,'figure_path',outputPath,'csv_path',csvPath);
end
