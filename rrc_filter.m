function [h, group_delay] = rrc_filter(alpha, span_symbols, sps)
%RRC_FILTER Generate a unit-energy root-raised-cosine FIR (toolbox-free).
%   [h, group_delay] = rrc_filter(alpha, span_symbols, sps)
%   alpha        - roll-off factor 0..1.
%   span_symbols - even FIR length in symbol periods.
%   sps          - integer samples per symbol.
%   h            - column vector with sum(abs(h).^2) = 1.
%   group_delay  - integer (length(h)-1)/2 in SAMPLES.
%
% Time t is normalized by the symbol period T (T=1). The special cases
% t=0 and abs(t)=1/(4*alpha) are evaluated using analytic limits. This
% avoids NaNs and Inf without relying on rcosdesign or any toolbox.

if ~(isscalar(alpha) && isfinite(alpha) && alpha >= 0 && alpha <= 1)
    error('rrc_filter:alpha', 'alpha must lie in [0,1].');
end
if ~(isscalar(span_symbols) && span_symbols > 0 && ...
     mod(span_symbols, 2) == 0 && span_symbols == floor(span_symbols))
    error('rrc_filter:span', 'span_symbols must be a positive even integer.');
end
if ~(isscalar(sps) && sps >= 2 && sps == floor(sps))
    error('rrc_filter:sps', 'sps must be an integer >= 2.');
end

t = (-span_symbols*sps/2 : span_symbols*sps/2).' / sps;
h = zeros(size(t));
eps_t = 1e-12;
if alpha == 0
    h = ones(size(t));
    nz = abs(t) > eps_t;
    h(nz) = sin(pi*t(nz)) ./ (pi*t(nz)); % sinc(t), implemented locally.
else
    at_zero = abs(t) < eps_t;
    at_sing = abs(abs(t) - 1/(4*alpha)) < eps_t;
    regular = ~(at_zero | at_sing);
    tr = t(regular);
    h(regular) = (sin(pi*tr*(1-alpha)) + ...
                  4*alpha*tr.*cos(pi*tr*(1+alpha))) ./ ...
                 (pi*tr.*(1-(4*alpha*tr).^2));
    h(at_zero) = 1 - alpha + 4*alpha/pi;
    h(at_sing) = alpha/sqrt(2) * ( ...
        (1+2/pi)*sin(pi/(4*alpha)) + ...
        (1-2/pi)*cos(pi/(4*alpha)) );
end
h = h / sqrt(sum(abs(h).^2));
group_delay = (numel(h)-1)/2;
end
