clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIRiT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load DATA
% load phantom.mat
load brain_8ch

%% reconstruction parameters	
kSize = [5,5];  % SPIRiT kernel size
nIterCG = 20; % number of iteration; phantom requires twice as much as the brain.
mask_type = 'unif'; % options are: 'unif','random4','random3'
CalibTyk = 0.01;  % Tykhonov regularization in the calibration
ReconTyk = 1e-5;  % Tykhovon regularization in the reconstruction (SPIRiT only)

im = ifft2c(DATA);

%% mask decision
switch mask_type
    case 'unif'
           mask = mask_unif;
    case 'random3'
            mask = mask_randm_x3;            
    case 'random4'
            mask = mask_randm_x4;
end

%% initilisation
[pe,fe,coils] = size(DATA);
[CalibSize, dcomp] = getCalibSize(mask);  % get size of calibration area(full sampling center size) from mask

DATA = DATA.*mask; % masking
DATAcomp = DATA.*dcomp;
scale_fctr = norm(DATAcomp(:))/sqrt(coils)/20;  % scale_fctr = 1; // idk why it exists
DATA = DATA/scale_fctr;
DATAcomp = DATAcomp/scale_fctr;

im_dc = ifft2c(DATAcomp); % undersampled image
im = im/scale_fctr; % reference image

%% SPIRiT
disp('SPIRiT reconstruction')
kCalib = crop(DATA,[CalibSize,coils]); % crop the data which is going to used for calibration
kernel = zeros([kSize,coils,coils]); % // kernel size for recon : 3-D (1 kernel recon 1 data point in one coil) => needs 1 more dimension for referring recon data position(each coil idx)
[AtA,] = corrMatrix(kCalib,kSize);
for n=1:coils
	kernel(:,:,:,n) = calibrate(AtA,kSize,coils,n,CalibTyk);
end
GOP = SPIRiT(kernel, 'fft',[fe,pe]);

% CG
x_init = zeros(size(DATA));
tic;
[res_cg,rmse] = cgsp(DATA,GOP,nIterCG,ReconTyk,x_init,'cg'); % cg(malfunction) or lsqr(work)
toc;

% [res_cg, RESVEC] = cgSPIRiT(DATA,GOP,nIterCG,ReconTyk, DATA);

%% termination
im_cgspirit = ifft2c(res_cg); im_cgspirit_err = im_cgspirit - im;

%% Display
im_cgspirit_sqr = sos(im_cgspirit); im_cgspirit_err_sqr = sos(im_cgspirit_err);

im_dc_sqr = sos(im_dc);
im_sqr = sos(im);

figure, imshow(cat(2,im_sqr,im_dc_sqr,im_cgspirit_sqr),[]);
title ('REFERENCE                           Zero-padded                               SPIRiT');
figure, imshow(im_cgspirit_err_sqr,[]); title ('Difference image');
figure, plot(rmse); xlabel('iteration'); ylabel('rmse');



