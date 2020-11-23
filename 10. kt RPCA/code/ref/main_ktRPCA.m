% Simple demo and implementation of the following paper:
% Tremoulheac et al, Dynamic MR image reconstruction--separation from 
% undersampled (k,t)-space via low-rank plus sparse prior, TMI 2014
%
% Benjamin Tremoulheac
% b.tremoulheac@cs.ucl.ac.uk
% July 2014

clear all;
close all;
clc;

load phan.mat;
load mask;

[nx, ny, nt] = size(X0);

%% Define operators

E = @(Z) FSSO(Z,mask,1);
Eh = @(Z) FSSO(Z,mask,2);
EhE = @(Z) Eh(E(Z));

T = @(Z) tempFT(Z,[nx ny nt],1);
Th = @(Z) tempFT(Z,[nx ny nt],2);

%% Measurements, ZF-IDFT reconstruction

y = E(X);
recidft = reshape(Eh(y),[nx ny nt]);

%% k-t RPCA parameters

% Regularisation parameter
mu = 200;

% Decomposition parameter
rho = 1;
lambda = rho*(max(nx*ny,nt))^(-1/2); 


%% Call k-t RPCA

optsktrpca.quiet = false;
optsktrpca.tol = 10e-5;
optsktrpca.maxIter = 100;
optsktrpca.datasize = [nx ny nt];
optsktrpca.E = E;
optsktrpca.Eh = Eh;
optsktrpca.EhE = EhE;
optsktrpca.T = T;
optsktrpca.Tt = Th;
optsktrpca.mu1 = mu;
optsktrpca.mu2 = mu*lambda;
outktrpca = ktRPCA_ADMM(y,optsktrpca);
recktrpca = reshape(outktrpca.recon,[nx ny nt]);
recktrpca_L = reshape(outktrpca.L,[nx ny nt]);
recktrpca_S = reshape(outktrpca.S,[nx ny nt]);


%% Error

fprintf('\nUnder-sampling ratio = %g (acc.factor = %g)\n',...
    nnz(mask)/numel(mask),numel(mask)/nnz(mask));

erridft = (norm(reshape(recidft,[nx*ny nt])-X0c,'fro')^2)/(norm(X0c,'fro')^2);
errktrpca = (norm(reshape(recktrpca,[nx*ny nt])-X0c,'fro')^2)/(norm(X0c,'fro')^2);
fprintf('ZF-IDFT = %g dB\n',-10*log10(erridft));
fprintf('k-t RPCA = %g dB\n',-10*log10(errktrpca));

%% Display

while 1 == 1
    for k = 1:nt
    subplot(151); imagesc(abs(X0(:,:,k)),[0 255]); title(['Original noiseless - ',num2str(k)]); axis image off;
    subplot(152); imagesc(abs(recidft(:,:,k)),[0 255]); title('ZF-IDFT'); axis image off;   
    subplot(153); imagesc(abs(recktrpca(:,:,k)),[0 255]); title('L+S'); axis image off;
    subplot(154); imagesc(abs(recktrpca_L(:,:,k)),[0 255]); title('L'); axis image off;
    subplot(155); imagesc(abs(recktrpca_S(:,:,k)),[0 255]);  title('S');axis image off;
    pause(0.1);
    clf;
    end
end
