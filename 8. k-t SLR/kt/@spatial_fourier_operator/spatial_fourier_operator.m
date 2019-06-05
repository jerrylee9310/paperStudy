function res = spatial_fourier_operator(mask)

res.adjoint = 0; % transpose flag
res.mask = mask;

res = class(res,'spatial_fourier_operator');