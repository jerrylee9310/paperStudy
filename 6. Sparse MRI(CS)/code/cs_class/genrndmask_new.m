function [mask_opt,pdf] = genrndmask_new(kdata, pwrfct, pctg)
[sx,sy] = size(kdata);
%% pdf generation
% dc point
r = floor(size(kdata)/2) + 1; 

% distance from dc point
[x,y] = meshgrid(1-r:sx-r,sy-r+1:-1:2-r);
d = x + 1i*y;
d = abs(d); d = d/max(d(:));

% probablility distribution matrix
pdf = (1-d).^pwrfct;

% condition check
pwr_samples = floor(sum(pdf(:))); % controllable (can be bigger)
pctg_samples = floor(pctg*sx*sy); % uncontrollable (fixed)
if pctg_samples < pwr_samples
    error('Infeasible case, increase powerfactor or percentage')
end

% bisect(pdf optimization)
maxlim = pctg_samples/(sx*sy);
minlim = 0;
while(1)
    offset = maxlim/2 + minlim/2;
    pdf_tmp = pdf + offset;
    pdf_tmp(pdf_tmp>1) = 1;
    expected_samples = floor(sum(pdf_tmp(:)));
    if expected_samples == pctg_samples
        pdf = pdf_tmp;
        break;
    elseif expected_samples > pctg_samples
        maxlim = offset;
    elseif expected_samples < pctg_samples
        minlim = offset;
    end
end

%% monte-carlo algorithm
monte_maxiter = 100;
rmse_opt = sx*sy;

for i = 1:monte_maxiter
    mask = rand(size(pdf))<pdf;
    interf = ifft2c(mask./pdf);
    interf_E = abs(interf(:));
    interf_rmse = sqrt(sum(interf_E) - max(interf_E));
    
    if rmse_opt > interf_rmse
        rmse_opt = interf_rmse;
        mask_opt = mask;
    end
end





