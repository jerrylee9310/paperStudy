clear all; close all; clc; clear classes;
%%
load aperiodic_pincat.mat
img = new;
kdata = fft2(img); 
kdata = kdata/max(abs(kdata(:))); % normalise

maskparams.type = 'r'; % radial
maskparams.nspoke = 24;
% maskparams.type = 'c'; % cartesian
% maskparams.lowfreq = 30;

mask = genmask(size(img),maskparams);
mask = fftshift(fftshift(mask,1),2);

%% parameters
% regularization weight
params.reg.nuclear = 3e-2;
params.reg.sparsity = 1e-3;

% penalty weight
params.constr.casorati = 3e-4;
params.constr.sparsity = 5e-8;

params.al.incfact.casorati = 15;
params.al.incfact.sparsity = 25;
params.al.maxitin = 10;
params.al.maxit = 10;
params.al.tol = 1e-100;

% types of sparsity-promote prior
% params.prior.type = 'W'; % wavelet
params.prior.type = 'D'; % TV
params.prior.numofprior = 3;
% params.prior.type = 'WD'; % wavelet + TV

%% function handlers
Fspat = spatial_fourier_operator(mask);
Spar = sparsity_prior_operator(params);

%% AL famework
ud = mask.*kdata;

[recon, cost] = alslr(Fspat,Spar,ud,params);

%% plotting
x_init = Fspat'*ud;
m = fftshift(fftshift(mask,1),2);

figure(1); colormap(gray);
subplot(3,3,1); imagesc(abs(img(:,:,24))); title('Gold standard'); 
subplot(3,3,2); imagesc(abs(x_init(:,:,24))); title('Undersampled data'); 
subplot(3,3,3); imagesc(abs(recon(:,:,24))); title('k-t SLR'); 
subplot(3,3,4); imagesc(abs(squeeze(img(60,:,:)))); title('Ground truth, x-t domain'); 
subplot(3,3,5); imagesc(abs(squeeze(x_init(60,:,:)))); title('Undersampled data, x-t domain'); 
subplot(3,3,6); imagesc(abs(squeeze(recon(60,:,:)))); title('k-t SLR, x-t domain'); 
subplot(3,3,7); imagesc(abs(m(:,:,1)));
subplot(3,3,8); imagesc(abs(squeeze(m(68,:,:))));
subplot(3,3,9); plot(cost,'-x'); title('cost function value'); xlabel('iteration');


















