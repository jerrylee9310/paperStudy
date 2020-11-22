clear all; close all; clc;
addpath(strcat(pwd,'/utils')); 
%% 1. load data
% load brainsos
load im1
img = im1;

% normalising
img = img./max(abs(img(:)));
kdata = fft2c(img);

% sampling mask
pwrfct = 3;
pctg = 0.2;
[mask,pdf] = genrndmask_new(kdata,pwrfct,pctg);

% initial image
kdata = kdata.*mask;
img_init = ifft2c((kdata./pdf));

%% 2.Parameters setting
% iteration parameter
it_para.maxiter = 20;
it_para.gradtol = 1e-30;
it_para.alpha = 0.05;
it_para.beta = 0.6;
it_para.tmax = 20;

% wavelet parameter
w_para.filter = 'db1';
w_para.level = 4;
w_para.smoothingfct = 1e-6; % [1e-15, 1e-6]
w_para.weight = 1e-3;
[nothing,w_para.idx] = wav2c(img,w_para.level,w_para.filter);
W = wavletoperator(img,w_para.filter,w_para.level);
w_para.operator = W;

% fourier parameter
FT = ftoperator(mask);
f_para.operator = FT;

% TV parameter
TV = tvoperator();
tv_para.smoothingfct = 1e-6; % [1e-15, 1e-6]
tv_para.weight = 1e-3;
tv_para.operator = TV;

%% 3. Nonlinear CG with line backtracking
tic
[rec,rmse] = nlcg_new(kdata,it_para, w_para, f_para, tv_para);
toc

%% 4. plotting
tmp = cat(2,img,abs(img_init),abs(rec));
% reconstruction image
figure(), imshow(abs(tmp),[]); title('REF, under, recon')
% % RMSE
% figure(), plot(rmse(2:end),'-x'); title('RMSE'); xlabel('iteration'); ylabel('rmse')
% % center line (phantom)
% figure(), plot(img(256,:)); hold on; plot(abs(img_init(256,:))); plot(abs(rec(256,:))); legend('reference', 'z-pad','recon');

%%
rmpath(strcat(pwd,'/utils')); 