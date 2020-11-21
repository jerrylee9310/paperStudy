function res = ifft2c(x)

% res = ifft2c(x)
% 
% orthonormal inverse 2D FFT
%
% (c) Michael Lustig 2005

res = 1/sqrt(length(x(:)))*fftshift(ifft2(ifftshift(x)));