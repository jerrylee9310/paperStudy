clear all; close all; clc;
psize = 256;
calisize = 20;
addNoiseSTD = 0.03; % add noise. Use only for phantom. Brain has enough noise!

[x,y] = meshgrid(linspace(0,1,psize));
% Generate fake Sensitivity maps
sMaps = cat(3,x.^2,1-x.^2,y.^2,1-y.^2);
% generate 4 coil phantom
imgs = repmat(phantom(psize),[1,1,4]).*sMaps;
DATA = fft2c(imgs);

% add noise
DATA = DATA+randn(size(DATA))*addNoiseSTD + 1i*randn(size(DATA))*addNoiseSTD;

% crop 20x20 window from the center of k-space for calibration
kCalib = crop(DATA,[calisize,calisize,4]);

%calibrate a kernel
kSize = [5,5];
coils = 4;
kernel = zeros([kSize,coils,coils]);
[AtA,] = corrMatrix(kCalib,kSize);
for n=1:coils
  kernel(:,:,:,n) = calibrate(AtA,kSize,coils,n,0.01);
end

% undersample by a factor of 2
mask = ones(size(DATA));
mask(1:2:end,2:2:end) = 0;
mask(2:2:end,1:2:end) = 0;
DATA = DATA.*mask;

%% iteration

% kernel operator
GOP = SPIRiT(kernel, 'fft',[psize,psize]);
% iterative method
x_init = zeros(size(DATA));
[res,rmse] = cgsp(DATA,GOP, 20, 1e-5, x_init,'cg'); % 'cg' or 'lsqr'

%%
figure, imshow(cat(2,sos(imgs), 2*sos(ifft2c(DATA)), sos(ifft2c(res))),[]);
title('full,  zero-fill,   recon')
figure, plot(rmse);