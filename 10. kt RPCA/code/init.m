% clear all; close all; clc;

load phantom.mat;
[nx,ny,nt] = size(img);

% maskparams.type = 'r'; % radial
% maskparams.nspoke = 6;
maskparams.type = 'c'; % cartesian
maskparams.lowfreq = 6;
maskparams.rfac = 0.08;

mask = genmask(size(img),maskparams);
mask = fftshift(fftshift(mask,1),2);

%% operators
E = @(x) FTS(x,mask,'forward');
Eh = @(x) FTS(x,mask,'adjoint');
EhE = @(x) Eh(E(x));

Ft = @(x) FTT(x,[nx ny nt],'forward');
Fth = @(x) FTT(x,[nx ny nt],'adjoint');

%% Measurements, ZF-IDFT reconstruction
d = E(img);
uimg = reshape(Eh(d),[nx ny nt]);

%% parameters
rho = 1; 
params.lambda = rho*(max(ny*nx,nt))^(-1/2);  % Decomposition parameter
params.mu = 200; % Regularisation parameter
params.penalty_LR = 1;
params.penalty_SP = 1;

params.tol = 1e-4;
params.maxit = 100;
params.cgtol = 1e-4;
params.cgit = 100;

params.dsize = [nx,ny,nt];
params.E = E;
params.Eh = Eh;
params.EhE = EhE;
params.Ft = Ft;
params.Fth = Fth;

%%
recon = ktRPCA(d,params);

%% plotting
% m = fftshift(fftshift(mask,1),2);
% 
% figure(); colormap(gray);
% subplot(3,3,1); imagesc(abs(img(:,:,24))); title('Gold standard'); 
% subplot(3,3,2); imagesc(abs(uimg(:,:,24))); title('Undersampled data'); 
% subplot(3,3,3); imagesc(abs(recon(:,:,24))); title('kt RPCA'); 
% subplot(3,3,4); imagesc(abs(squeeze(img(60,:,:)))); title('Ground truth, x-t domain'); 
% subplot(3,3,5); imagesc(abs(squeeze(uimg(60,:,:)))); title('Undersampled data, x-t domain'); 
% subplot(3,3,6); imagesc(abs(squeeze(recon(60,:,:)))); title('k-t PS, x-t domain'); 
% subplot(3,3,7); imagesc(abs(m(:,:,1))); xlabel('Readout'); ylabel('Phase Encoding');
% subplot(3,3,8); imagesc(abs(squeeze(m(:,60,:)))); xlabel('Time frame'); ylabel('Phase Encoding');
fprintf('Sampling ratio : %f\n', sum(mask(:))/numel(mask));

m = fftshift(fftshift(mask,1),2);
iimg = img/max(abs(img(:))); uuimg = uimg/max(abs(uimg(:)));
rimg = recon/max(abs(recon(:)));
% iimg = abs(img/max(abs(img(:)))); uuimg = abs(uimg/max(abs(uimg(:))));
% rimg = abs(recon/max(abs(recon(:))));
uerr = abs(iimg - uuimg);
rerr = abs(iimg - rimg);

tm = 24;
slc = 60;

goldstd = [iimg(:,:,tm); angle(iimg(:,:,tm))/3; imresize(squeeze(iimg(slc,:,:)),[128 128])];
initd = [uuimg(:,:,tm); angle(uuimg(:,:,tm))/3; imresize(squeeze(uuimg(slc,:,:)),[128 128])];
recd = [rimg(:,:,tm); angle(rimg(:,:,tm))/3; imresize(squeeze(rimg(slc,:,:)),[128 128])];
mmask =  [m(:,:,1) imresize(squeeze(m(:,60,:)),[128 128])];

figure(); imshow(cat(2,goldstd,initd,recd),[0 0.8]);
figure(); imshow(mmask,[]);





