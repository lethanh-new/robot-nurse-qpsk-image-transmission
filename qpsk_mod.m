function s = qpsk_mod(bits)
%QPSK_MOD Gray-labeled unit-energy QPSK modulator (complex baseband).
% Input : bits - even-length vector of binary values (0 or 1).
% Output: s    - complex QPSK symbols, one sample per symbol, Es = 1.
%
% Mapping, first bit controls in-phase (I), second quadrature (Q):
%   00 -> (+1 + j)/sqrt(2)      01 -> (+1 - j)/sqrt(2)
%   10 -> (-1 + j)/sqrt(2)      11 -> (-1 - j)/sqrt(2)
% Adjacent points on the constellation differ by just one bit (Gray).

if ~isvector(bits) || mod(numel(bits), 2) ~= 0
    error('qpsk_mod:InputLength', 'bits must be an even-length vector.');
end
if ~isnumeric(bits) && ~islogical(bits)
    error('qpsk_mod:InputType', 'bits must be numeric or logical.');
end
if any(bits(:) ~= 0 & bits(:) ~= 1)
    error('qpsk_mod:InputValues', 'bits must contain only 0 or 1.');
end

pairs = reshape(double(bits(:)), 2, []);  % Row 1: I bit; row 2: Q bit
s = ((1 - 2*pairs(1,:)) + 1j*(1 - 2*pairs(2,:))).' / sqrt(2);
end
