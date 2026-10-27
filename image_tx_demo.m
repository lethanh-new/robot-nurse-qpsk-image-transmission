function results = image_tx_demo(image_path)
%IMAGE_TX_DEMO Transmit a 128x128 grayscale image through the full QPSK RF link.
%   results = image_tx_demo()
%   results = image_tx_demo('images/06_ecg_flatline_128.png')
%
% Chain: image -> bitstream -> Gray QPSK -> upsample/RRC -> real RF at 40 kHz
%      -> AWGN -> coherent downconversion -> matched RRC -> hard decision
%      -> bitstream -> reconstructed image.
%
% Outputs:
%   - figures/image_reconstruction_grid.png
%   - figures/image_compare_XXdB.png (one figure per Eb/N0 point)
%   - figures/image_restored_XXdB.png (raw reconstructed image)
%   - figures/image_results_table.csv
%
% The FFT/Fourier transform was used in earlier waveform diagnostics to show
% spectral occupancy and to implement long FIR filtering efficiently. Here we
% focus on end-to-end image transport quality metrics: BER, MSE, and PSNR.

p = params();
rng(p.rng_seed, 'twister');

% Resolve paths relative to this source file, not MATLAB's current folder.
% An absolute output directory also avoids ambiguity from similarly named
% directories elsewhere on the MATLAB search path.
projectDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(projectDir, p.figure_dir);
if nargin < 1 || isempty(image_path)
    image_path = fullfile(projectDir, p.image_input_path);
elseif exist(image_path, 'file') ~= 2
    % Support relative custom paths even when MATLAB starts elsewhere.
    image_path = fullfile(projectDir, char(image_path));
end

% Test PNG writing BEFORE running the expensive RF simulation.
if exist(outputDir, 'dir') ~= 7
    [made, why] = mkdir(outputDir);
    if ~made && exist(outputDir, 'dir') ~= 7
        error('image_tx_demo:CreateDir', ...
            'Could not create output directory: %s\nReason: %s', outputDir, why);
    end
end
probeFile = [tempname(outputDir), '.png'];
try
    imwrite(uint8(zeros(2,2)), probeFile);
    delete(probeFile);
catch ME
    error('image_tx_demo:NotWritable', ...
        ['Cannot write PNG files in: %s\n' ...
         'MATLAB reported: %s\n' ...
         'Check folder permissions or change p.figure_dir in params.m.'], ...
         outputDir, ME.message);
end
fprintf('Verified PNG output folder: %s\n', outputDir);

[txBits, imgTx] = image_to_bits(image_path);
txSymbols = qpsk_mod(txBits);
EbN0_list = p.image_EbN0_dB(:).';
nPts = numel(EbN0_list);

numBits = zeros(nPts,1);
numErrors = zeros(nPts,1);
BER_img = zeros(nPts,1);
MSE_img = zeros(nPts,1);
PSNR_img = zeros(nPts,1);
restoredFiles = strings(nPts,1);
compareFiles = strings(nPts,1);
imgRx_all = cell(nPts,1);

[~, baseName, ~] = fileparts(image_path);
fprintf('\nROBOT NURSE -> MONITORING STATION | IMAGE TRANSMISSION\n');
fprintf('Input image: %s\n', image_path);
fprintf('Image size : %d x %d pixels, %d bits\n', size(imgTx,1), size(imgTx,2), numel(txBits));
fprintf('Mode       : Gray QPSK, RRC shaped, fc = %.0f Hz, fs = %.0f Hz, sps = %d\n', ...
    p.fc_Hz, p.fs_Hz, p.sps);
fprintf('rng(%d), Eb/N0 points = [%s] dB\n\n', p.rng_seed, sprintf('%g ', EbN0_list));

for k = 1:nPts
    EbN0_dB = EbN0_list(k);
    rxSymbols = qpsk_passband_link(txSymbols, EbN0_dB, p, true);
    rxBits = qpsk_demod(rxSymbols);
    rxBits = rxBits(1:numel(txBits));
    imgRx = bits_to_image(rxBits, size(imgTx));
    imgRx_all{k} = imgRx;

    numBits(k) = numel(txBits);
    numErrors(k) = sum(rxBits ~= txBits);
    BER_img(k) = numErrors(k) / numBits(k);
    [MSE_img(k), PSNR_img(k)] = image_metrics(imgTx, imgRx);

    restoredName = sprintf('%s%s_%ddB.png', p.image_restored_prefix, baseName, EbN0_dB);
    restoredPath = fullfile(outputDir, restoredName);
    imwrite(imgRx, restoredPath);
    restoredFiles(k) = string(restoredPath);

    compareName = sprintf('%s%s_%ddB.png', p.image_compare_prefix, baseName, EbN0_dB);
    comparePath = fullfile(outputDir, compareName);
    save_comparison_figure(imgTx, imgRx, EbN0_dB, BER_img(k), MSE_img(k), PSNR_img(k), comparePath, p.png_dpi);
    compareFiles(k) = string(comparePath);

    fprintf('Processed Eb/N0 = %2d dB: bit errors = %5d, BER = %.6e, MSE = %.3f, PSNR = %.3f dB\n', ...
        EbN0_dB, numErrors(k), BER_img(k), MSE_img(k), PSNR_img(k));
end

results = table(EbN0_list(:), numBits, numErrors, BER_img, MSE_img, PSNR_img, restoredFiles, compareFiles, ...
    'VariableNames', {'EbN0_dB','num_bits','num_errors','BER_image','MSE','PSNR_dB','restored_png','compare_png'});

fprintf('\nEbN0(dB) |    Bits | Errors |    BER_image |       MSE |  PSNR(dB)\n');
fprintf('-------------------------------------------------------------------\n');
for k = 1:nPts
    fprintf('%8.1f | %7d | %6d | %12.6e | %9.3f | %9.3f\n', ...
        results.EbN0_dB(k), results.num_bits(k), results.num_errors(k), ...
        results.BER_image(k), results.MSE(k), results.PSNR_dB(k));
end

csvPath = fullfile(outputDir, p.image_table_csv);
writetable(results, csvPath);
fprintf('Saved: %s\n', csvPath);

% Summary grid for report: original + one reconstructed image per Eb/N0.
gridPath = fullfile(outputDir, p.image_grid_png);
save_grid_figure(imgTx, imgRx_all, EbN0_list, BER_img, PSNR_img, gridPath, p.png_dpi);
fprintf('Saved: %s\n', gridPath);
end

function save_comparison_figure(imgTx, imgRx, EbN0_dB, berVal, mseVal, psnrVal, outPath, dpi)
fig = figure('Visible', 'off', 'Color', 'w', 'Name', 'Image comparison');
subplot(1,2,1);
imagesc(imgTx, [0 255]);
axis image off;
colormap(gray(256));
title('Original 128x128 grayscale');

subplot(1,2,2);
imagesc(imgRx, [0 255]);
axis image off;
colormap(gray(256));
title(sprintf('Reconstructed @ E_b/N_0 = %d dB\nBER = %.3e, MSE = %.2f, PSNR = %.2f dB', ...
    EbN0_dB, berVal, mseVal, psnrVal));

sgtitle(sprintf('Robot Nurse image transmission: original vs reconstructed (%d dB)', EbN0_dB));
print(fig, outPath, '-dpng', sprintf('-r%d', dpi));
close(fig);
end

function save_grid_figure(imgTx, imgRx_all, EbN0_list, berList, psnrList, outPath, dpi)
fig = figure('Visible', 'off', 'Color', 'w', 'Name', 'Image reconstruction grid');
numPanels = numel(EbN0_list) + 1;
rows = 2;
cols = ceil(numPanels / rows);

subplot(rows, cols, 1);
imagesc(imgTx, [0 255]);
axis image off;
colormap(gray(256));
title('Original image');

for ii = 1:numel(EbN0_list)
    subplot(rows, cols, ii+1);
    imagesc(imgRx_all{ii}, [0 255]);
    axis image off;
    colormap(gray(256));
    title(sprintf('%d dB\nBER = %.2e\nPSNR = %.2f dB', EbN0_list(ii), berList(ii), psnrList(ii)));
end

sgtitle('Robot Nurse -> Monitoring Station: reconstructed images through QPSK link');
print(fig, outPath, '-dpng', sprintf('-r%d', dpi));
close(fig);
end
