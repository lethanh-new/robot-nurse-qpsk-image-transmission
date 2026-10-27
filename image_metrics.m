function [mse_val, psnr_val] = image_metrics(img_ref, img_test)
%IMAGE_METRICS Compute MSE and PSNR for uint8 grayscale images.

if ~isequal(size(img_ref), size(img_test))
    error('image_metrics:Size', 'Input images must have the same size.');
end
err = double(img_ref) - double(img_test);
mse_val = mean(err(:).^2);
if mse_val == 0
    psnr_val = Inf;
else
    psnr_val = 10 * log10((255^2) / mse_val);
end
end
