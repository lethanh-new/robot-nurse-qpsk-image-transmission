function results = image_tx()
%IMAGE_TX Synthetic 64x64 PCM/JPEG image transmission through full QPSK RF.
%   results = image_tx();
%
% Only MATLAB base functions are required (no comms or image toolboxes).
% Reuses GĐ2 pipeline qpsk_mod, qpsk_passband_link, qpsk_demod:
%  bits -> Gray QPSK (Es=1) -> x16 upsample -> Tx RRC -> real 40 kHz RF
%       -> real AWGN -> coherent downconversion -> matched Rx RRC
%       -> sample after Tx+Rx group delay -> bits -> image / JPEG bytes.
%
% Encoding comparisons:
%   PCM8: 64*64*8 raw bits (MSB first), directly rebuild uint8 pixels.
%   JPEG: actual compressed file bytes (headers AND entropy data) are sent;
%         the receiver writes received bytes back to a JPEG and attempts
%         imread. Damaged headers or coded data may make decoding fail.
%         Decode failure is recorded, NEVER silently replaced/estimated.
%
% MSE = mean((double(original)-double(reconstruction)).^2).
% PSNR = 10*log10(255^2/MSE), Inf only if MSE=0.
% File duration = payload bits / Rb; excludes RF guard/frame overhead.
% FFT/Fourier transform in GĐ2 implements zero-padded RRC FIR convolution
% and diagnoses occupied bandwidth (waveform_demo.m); this script compares
% reconstruction, not the frequency spectrum.
%
% Fixed settings (image size, SNR, JPEG quality, output names) in params.m.
% Reproducible channel noise uses rng(1) from params.m.

p = params();
rng(p.rng_seed, 'twister');
projectDir = fileparts(mfilename('fullpath'));
outDir = fullfile(projectDir, p.figure_dir);
if exist(outDir, 'dir') ~= 7
    [ok, msg] = mkdir(outDir);
    if ~ok && exist(outDir, 'dir') ~= 7
        error('image_tx:OutputDir', 'Cannot create %s: %s', outDir, msg);
    end
end

% Verify we can write BEFORE performing waveform simulations.
probePath = [tempname(outDir), '.png'];
try
    imwrite(uint8(zeros(2,2)), probePath);
    delete(probePath);
catch ME
    error('image_tx:WritePermission', 'Cannot write into %s: %s', outDir, ME.message);
end

imgOriginal = make_synthetic_image64(p);
imageShape = size(imgOriginal);
if ~isa(imgOriginal,'uint8') || ~isequal(imageShape, [64 64])
    error('image_tx:SyntheticInput', 'Expected synthetic uint8 grayscale 64x64.');
end
imwrite(imgOriginal, fullfile(outDir, p.image64_png));

% Encode PCM8 pixels in MATLAB column-major order, bit significance MSB-first.
pcmBytes = imgOriginal(:);
pcmBits = bytes_to_bits(pcmBytes);
if ~isequal(reshape(bits_to_bytes(pcmBits), imageShape), imgOriginal)
    error('image_tx:PCMEncode', 'PCM8 lossless byte round-trip failed.');
end

% Encode the synthetic image as a REAL lossy JPEG file (not a fake bit count).
sourceJPEG = fullfile(outDir, p.image64_jpeg);
imwrite(imgOriginal, sourceJPEG, 'jpg', 'Quality', p.jpeg_quality);
readFid = fopen(sourceJPEG, 'rb');
if readFid < 0
    error('image_tx:JPEGRead', 'Cannot read compressed source: %s', sourceJPEG);
end
jpegBytes = fread(readFid, Inf, '*uint8');
fclose(readFid);
jpegBits = bytes_to_bits(jpegBytes);

% Decode error-free JPEG to quantify intrinsic lossy compression.
baseJPEG = imread(sourceJPEG);
if ndims(baseJPEG) ~= 2 || ~isequal(size(baseJPEG), imageShape)
    error('image_tx:JPEGDimensions', 'JPEG baseline is not 64x64 grayscale.');
end
[baseMSE, basePSNR] = image_metrics(imgOriginal, baseJPEG);

Rb = (p.fs_Hz/p.sps) * p.bits_per_sym;
EbN0dB = p.image64_ebn0_dB(:).';
numSnr = numel(EbN0dB);
modeNames = {'PCM8', 'JPEG'};
modeBits = {pcmBits, jpegBits};
numModes = numel(modeNames);
nRows = numSnr * numModes;

modeOut = strings(nRows,1);
snrOut = zeros(nRows,1);
bitsOut = zeros(nRows,1);
errorsOut = zeros(nRows,1);
berOut = zeros(nRows,1);
mseOut = nan(nRows,1);
psnrOut = nan(nRows,1);
secsOut = zeros(nRows,1);
statusOut = strings(nRows,1);
sourceBytesOut = zeros(nRows,1);
imagePanels = cell(numModes, numSnr);

fprintf('\nROBOT NURSE | SYNTHETIC 64x64 PCM8 vs JPEG, FULL RRC-RF QPSK\n');
fprintf('Image: 64x64 uint8, synthetic monitor graphic (no patient images).\n');
fprintf('fs = %.0f Hz, fc = %.0f Hz, sps = %d, Rb = %.0f bit/s\n', ...
    p.fs_Hz, p.fc_Hz, p.sps, Rb);
fprintf('RRC alpha = %.2f, span = %d symbols, JPEG quality = %d\n', ...
    p.rrc_alpha, p.rrc_span_symbols, p.jpeg_quality);
fprintf('rng(%d), Eb/N0 = [%s] dB\n', p.rng_seed, sprintf('%g ',EbN0dB));
fprintf('PCM8 size = %d bits (%d bytes), payload time = %.6f s\n', ...
    numel(pcmBits), numel(pcmBytes), numel(pcmBits)/Rb);
fprintf('JPEG size = %d bits (%d bytes), payload time = %.6f s\n', ...
    numel(jpegBits), numel(jpegBytes), numel(jpegBits)/Rb);
fprintf('JPEG compression ratio vs PCM8 = %.3f\n', numel(pcmBits)/numel(jpegBits));
fprintf('Lossless channel + lossy JPEG baseline: MSE = %.3f, PSNR = %.3f dB\n', ...
    baseMSE, basePSNR);
fprintf('Payload time excludes fixed RRC guard symbols and any packet overhead.\n\n');

for m = 1:numModes
    txBits = modeBits{m};
    txSymbols = qpsk_mod(txBits);
    for k = 1:numSnr
        row = (m-1)*numSnr+k;
        thisSnr = EbN0dB(k);
        rxSymbols = qpsk_passband_link(txSymbols, thisSnr, p, true);
        rxBits = qpsk_demod(rxSymbols);
        if numel(rxBits) ~= numel(txBits)
            error('image_tx:Frame', 'Receive bitstream length does not match payload.');
        end
        nErr = sum(rxBits ~= txBits);
        rxBytes = bits_to_bytes(rxBits);

        modeOut(row) = string(modeNames{m});
        snrOut(row) = thisSnr;
        bitsOut(row) = numel(txBits);
        errorsOut(row) = nErr;
        berOut(row) = nErr / numel(txBits);
        secsOut(row) = numel(txBits)/Rb;
        sourceBytesOut(row) = numel(txBits)/8;

        decoded = [];
        if m == 1
            decoded = reshape(rxBytes, imageShape); % PCM8 never depends on a parser.
            statusOut(row) = "DECODE_OK";
        else
            % WARNING: JPEG bit errors can corrupt markers, code tables,
            % variable-length Huffman data or stream termination.
            receivedJPEG = fullfile(outDir, sprintf('%s%ddB.jpg', ...
                p.image64_rx_jpeg_prefix, thisSnr));
            try
                writeFid = fopen(receivedJPEG, 'wb');
                if writeFid < 0
                    error('image_tx:JPEGWrite', 'Could not create JPEG byte file.');
                end
                count = fwrite(writeFid, rxBytes, 'uint8');
                fclose(writeFid);
                if count ~= numel(rxBytes)
                    error('image_tx:JPEGIncomplete', 'JPEG bytes were not fully written.');
                end
                lastwarn('');
                decoded = imread(receivedJPEG);
                [warningText, ~] = lastwarn();
                if ndims(decoded) ~= 2 || ~isequal(size(decoded), imageShape)
                    error('image_tx:BadJPEG', 'Decoded JPEG is not 64x64 grayscale.');
                end
                if ~isa(decoded, 'uint8')
                    error('image_tx:BadJPEGType', 'JPEG decoded with unexpected pixel type.');
                end
                if isempty(warningText)
                    statusOut(row) = "DECODE_OK";
                else
                    statusOut(row) = "DECODE_WARN";
                end
            catch ME
                decoded = [];
                statusOut(row) = "DECODE_FAIL";
                fprintf('  JPEG @ %g dB could not be decoded: %s\n', thisSnr, ME.message);
            end
        end

        if ~isempty(decoded)
            [mseOut(row), psnrOut(row)] = image_metrics(imgOriginal, decoded);
            recoveredName = sprintf('%s%s_%ddB.png', ...
                p.image64_recovered_prefix, lower(modeNames{m}), thisSnr);
            imwrite(decoded, fullfile(outDir, recoveredName));
        end
        imagePanels{m,k} = decoded;
        fprintf('%-5s | %2g dB | %7d bits | %5d errors | BER %.6e | MSE %10.3f | PSNR %9.3f dB | %s\n', ...
            modeNames{m}, thisSnr, bitsOut(row), nErr, berOut(row), ...
            mseOut(row), psnrOut(row), char(statusOut(row)));
    end
end

results = table(modeOut, snrOut, bitsOut, errorsOut, berOut, mseOut, psnrOut, ...
    secsOut, sourceBytesOut, statusOut, ...
    'VariableNames', {'encoding','EbN0_dB','num_bits','num_errors','BER_image', ...
    'MSE','PSNR_dB','payload_time_s','source_bytes','decode_status'});

fprintf('\nCOMPLETE NUMERICAL TABLE (all 3 Eb/N0 levels x 2 encodings)\n');
disp(results);
outCSV = fullfile(outDir, p.image64_csv);
writetable(results, outCSV);
fprintf('Saved: %s\n', outCSV);

% Two 4-panel figures, same scale/intensity and labeled pixel axes.
for m = 1:numModes
    fig = figure('Color','w','Name',sprintf('%s: original vs received',modeNames{m}));
    for panel = 1:(numSnr+1)
        subplot(2,2,panel);
        if panel == 1
            imagesc(imgOriginal, [0 255]);
            title('Original synthetic 64x64');
        else
            k = panel-1;
            row = (m-1)*numSnr+k;
            if isempty(imagePanels{m,k})
                imagesc(zeros(imageShape,'uint8'), [0 255]);
                text(p.image64_side/2,p.image64_side/2, 'JPEG DECODE FAIL', ...
                    'Color','w','FontWeight','bold','HorizontalAlignment','center', ...
                    'FontSize',9);
                title(sprintf('%g dB | BER %.2e | DECODE FAIL', EbN0dB(k), berOut(row)));
            else
                imagesc(imagePanels{m,k}, [0 255]);
                title(sprintf('%g dB | BER %.3e | PSNR %.2f dB', ...
                    EbN0dB(k), berOut(row), psnrOut(row)));
            end
        end
        axis image;
        xlim([1 p.image64_side]); ylim([1 p.image64_side]);
        colormap(gray(256));
        xlabel('Column (pixels)'); ylabel('Row (pixels)');
    end
    if m == 1
        figTitle = sprintf('PCM8 raw image | %d bit/s | RRC QPSK RF', Rb);
        figFile = p.image64_pcm_fig;
    else
        figTitle = sprintf('JPEG (quality %d) | %d bit/s | RRC QPSK RF', p.jpeg_quality, Rb);
        figFile = p.image64_jpeg_fig;
    end
    sgtitle(figTitle);
    outPNG = fullfile(outDir, figFile);
    print(fig, outPNG, '-dpng', sprintf('-r%d',p.png_dpi));
    fprintf('Saved: %s\n', outPNG);
end
fprintf('\nWARNING: JPEG is inherently lossy. At BER=0 its PSNR can be finite.\n');
fprintf('At JPEG DECODE_FAIL, MSE/PSNR are NaN (not zero and not invented).\n');
end
