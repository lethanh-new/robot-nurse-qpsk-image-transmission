function ok = image_self_check()
%IMAGE_SELF_CHECK Verify 128x128 grayscale -> bits -> image round-trip.
% Also verifies the Gray-QPSK mapper and hard-demapper without channel noise.
% Requires no MATLAB toolbox. Prints one numerical test row per image.

p = params();
imageNames = dir(fullfile(p.image_dir, '*.png'));
if isempty(imageNames)
    error('image_self_check:NoImages', 'No PNG images in: %s', p.image_dir);
end

fprintf('IMAGE BITSTREAM SELF CHECK (128x128, uint8, MSB first)\n');
fprintf('%-35s | %8s | %10s | %12s\n', ...
    'Image', 'Bits', 'Pix errors', 'QPSK bit errors');
fprintf('%s\n', repmat('-', 1, 76));
ok = true;
for k = 1:numel(imageNames)
    name = imageNames(k).name;
    [bits, reference] = image_to_bits(fullfile(p.image_dir, name));
    reconstructed = bits_to_image(bits, size(reference));
    pixelErrors = sum(reference(:) ~= reconstructed(:));
    bitErrors = sum(qpsk_demod(qpsk_mod(bits)) ~= bits);
    fprintf('%-35s | %8d | %10d | %12d\n', ...
        name, numel(bits), pixelErrors, bitErrors);
    ok = ok && numel(bits) == p.image_bits && pixelErrors == 0 && bitErrors == 0;
end
if ok
    fprintf('PASS: all images round-trip exactly.\n');
else
    error('image_self_check:Failed', 'One or more round-trip checks failed.');
end
end
