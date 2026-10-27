function img = make_synthetic_image64(p)
%MAKE_SYNTHETIC_IMAGE64 Generate a synthetic 64x64 uint8 monitor pattern.
% No real patient images, identity data, or physiological recordings.
% The illustration contains a generated trace, geometric shapes and bars.
% All experiment settings, including image size, live in params.m.

N = p.image64_side;
if N ~= 64
    error('make_synthetic_image64:Size', 'This lab specifies 64x64 pixels.');
end
[x, y] = meshgrid(1:N, 1:N);
img = uint8(18 + round(12 * (x-1) / (N-1)) + 2*mod(floor((y-1)/8),2));

% Top monitor header and three synthetic indicator dots.
img(4:9, 4:61) = 66;
for cx = [13, 32, 51]
    mask = (x-cx).^2 + (y-6).^2 <= 4;
    img(mask) = uint8(215);
end

% A purely synthetic trace, NOT patient ECG or a clinical measurement.
img(14:37, 4:61) = 34;
for col = 5:60
    q = mod(col-5, 24);
    if q < 9
        row = 28;
    elseif q == 9
        row = 31;
    elseif q == 10
        row = 17;
    elseif q == 11
        row = 34;
    elseif q == 12
        row = 27;
    else
        row = 28 + round(2*sin((q-12)*pi/11));
    end
    img(max(14,row-1):min(37,row+1),col) = uint8(230);
end

% Synthetic strips and bars to make lossy JPEG artifacts visible.
img(42:44, 5:60) = 95;
for k = 1:5
    left = 7 + 11*(k-1);
    barheight = 3 + 2*k;
    img((59-barheight):59, left:(left+5)) = uint8(80 + 28*k);
end
img(61:62, 5:60) = uint8(120);
end
