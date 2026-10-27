function [bits, img] = image_to_bits(image_path)
%IMAGE_TO_BITS Read a grayscale 128x128 image and serialize to bits.
%   [bits, img] = image_to_bits(image_path)
% Output bits are column vector, MSB-first inside each uint8 pixel.
% This function uses only MATLAB base functionality.

if ~(ischar(image_path) || (isstring(image_path) && isscalar(image_path)))
    error('image_to_bits:Path', 'image_path must be a character vector or string scalar.');
end
image_path = char(image_path);
if exist(image_path, 'file') ~= 2
    error('image_to_bits:MissingFile', 'Image file not found: %s', image_path);
end

img = imread(image_path);
if ndims(img) == 3
    img = uint8(round(0.2989*double(img(:,:,1)) + 0.5870*double(img(:,:,2)) + 0.1140*double(img(:,:,3))));
end
if ~isa(img, 'uint8')
    img = uint8(img);
end
if ~isequal(size(img), [128 128])
    error('image_to_bits:Size', ...
        'Expected a 128x128 image. Use the prepared dataset or resize externally first.');
end

pixels = img(:);
% bitget needs a scalar bit-position or arrays of equal size on many
% MATLAB releases. Extract each bit plane separately (MSB to LSB).
bitMat = zeros(numel(pixels), 8, 'uint8');
for k = 1:8
    bitMat(:, k) = bitget(pixels, 9-k);
end
bits = reshape(bitMat.', [], 1);
end
