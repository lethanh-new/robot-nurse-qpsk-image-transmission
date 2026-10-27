function img = bits_to_image(bits, image_size)
%BITS_TO_IMAGE Reconstruct a uint8 grayscale image from an MSB-first bitstream.
%   img = bits_to_image(bits, [rows cols])

if nargin < 2
    image_size = [128 128];
end
if ~(isnumeric(image_size) && numel(image_size) == 2)
    error('bits_to_image:Size', 'image_size must be [rows cols].');
end
if mod(numel(bits), 8) ~= 0
    error('bits_to_image:Length', 'Number of bits must be a multiple of 8.');
end
if any(bits(:) ~= 0 & bits(:) ~= 1)
    error('bits_to_image:Values', 'bits must contain only 0 or 1.');
end

b = reshape(double(bits(:)), 8, []).';
weights = [128; 64; 32; 16; 8; 4; 2; 1];
pixels = uint8(b * weights);
if numel(pixels) ~= prod(image_size)
    error('bits_to_image:Count', 'Bit count does not match requested image size.');
end
img = reshape(pixels, image_size);
end
