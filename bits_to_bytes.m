function byteData = bits_to_bytes(bits)
%BITS_TO_BYTES Convert an MSB-first binary vector back into uint8 bytes.

if mod(numel(bits), 8) ~= 0
    error('bits_to_bytes:Length', 'Bitstream length must be divisible by 8.');
end
if any(bits(:) ~= 0 & bits(:) ~= 1)
    error('bits_to_bytes:Values', 'bits must contain only binary values.');
end
bits8 = reshape(double(bits(:)), 8, []).';
weights = [128; 64; 32; 16; 8; 4; 2; 1];
byteData = uint8(bits8 * weights);
end
