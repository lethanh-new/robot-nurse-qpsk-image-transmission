function bits = bytes_to_bits(byteData)
%BYTES_TO_BITS Serialize uint8 bytes MSB first to a binary column vector.
% Compatible with older MATLAB: bitget() is called with scalar bit index.

if ~isa(byteData, 'uint8')
    error('bytes_to_bits:Type', 'Input must be uint8.');
end
byteData = byteData(:);
bitPlanes = zeros(numel(byteData), 8, 'uint8');
for b = 1:8
    bitPlanes(:,b) = bitget(byteData, 9-b);
end
bits = reshape(bitPlanes.', [], 1);
end
