function bits_hat = qpsk_demod(r)
%QPSK_DEMOD Hard-decision Gray QPSK demodulator (complex baseband).
% Input : r        - complex received symbols (one sample/symbol).
% Output: bits_hat - column vector [Ibit1; Qbit1; Ibit2; Qbit2; ...].
% Ideal timing/carrier synchronization and zero phase offset are assumed.
% Decision boundaries are Re{r} = 0 and Im{r} = 0.

if ~isnumeric(r) || ~isvector(r)
    error('qpsk_demod:InputType', 'r must be a numeric vector.');
end
r = r(:);
bits_hat = zeros(2*numel(r), 1);
bits_hat(1:2:end) = real(r) < 0;
bits_hat(2:2:end) = imag(r) < 0;
end
