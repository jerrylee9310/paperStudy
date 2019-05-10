clear all; close all; clc;

load Brain2D_vertex;
DATA = DATA./max(abs(DATA(:))); % normalise on k-space
[nx,ny,nc] = size(DATA);
img = sos(ifft2(DATA));

% coil sensitivity
smap = getSmap(DATA);

% Data sampling
load SP
SP = fftshift(SP);
data = DATA.*SP;

% cf. indirect method
% load SP
% SP = fftshift(SP);
% data = fft2(img.*smap).*SP;

%% regularizor
% Type of sparsifying transform
% 	params.Operator = 'FD'; % finite differences
% 	params.Operator = 'W'; % Redundant wavelet transform
	params.Operator = 'WFD'; % Redundant wavelet transform & finite diff.

%% wavelet

params.Wavelet.wname = 'haar';
params.Wavelet.nlev = 2;
params.Wavelet.dnum = 3*params.Wavelet.nlev+1;

% Wavelet filters
[lod, hid, lor, hir] = wfilters(params.Wavelet.wname);

% Normalize filters so as to avoid a product with 0.5 during inverse undecimated wavelet transform 
params.Wavelet.lod = lod/sqrt(2); 
params.Wavelet.hid = hid/sqrt(2); 
params.Wavelet.lor = lor/sqrt(2); 
params.Wavelet.hir = hir/sqrt(2); 

%% parameters
params.nx = nx;
params.ny = ny;
params.nc = nc;
params.smap = smap;
params.mask = SP;

params.lambda.tv = 2e-8;
params.lambda.wavelet = 1.5e-7;

RR = compute_RR(params); % FFT of R^T * R
params.AL.RR = RR;

params.maxitr = 10;
params.maxitr_in = 10;

params.AL.mu = 5e-5;
params.AL.nu1 = 1;
params.AL.nu2 = 1;

% Inverse of some matrices required for solving sub-problems
iPpmu = 1 ./ (SP + params.AL.mu); % Freq. Response of (F'F + mu * I)^-1
iRpnu2nu1 = 1./(RR + (params.AL.nu2/params.AL.nu1)); % Freq. Response of (R'R + nu2/nu1 * I)^-1
iSpnu2 = 1./(sum(abs(smap).^2,3) + params.AL.nu2); % Inverse of (S'S + nu2 * I)
params.AL.iPpmu = iPpmu;
params.AL.iSpnu2= iSpnu2;
params.AL.iRpnu2nu1 = iRpnu2nu1;

eta20 = zeros(nx,ny,nc);
eta21 = zeros(nx,ny);
switch(params.Operator)
    case{'FD'} % Finite difference
        eta22 = zeros(nx,ny,2);
    case{'W'} % Wavelets
        eta22 = zeros(nx,ny,3*params.Wavelet.nlev+1);
    case{'WFD'}
        eta22 = zeros(nx,ny,2+3*params.Wavelet.nlev+1);
end

eta = cat(3,eta20,eta21,eta22);
u = zeros(size(eta));

%%
Fhd = ifft2(data);
x0 = sqrt(sum(abs(ifft2(data)).^2, 3));

[x_est, u0, u1, u2] = SENSERecon_ALP2(Fhd, eta, u, x0, params);
% x_est = SENSERecon_ALP2_CG(Fhd, eta, u, x0, params);

mshow(x_est);

