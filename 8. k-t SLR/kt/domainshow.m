clear all; clc; close all;

load('aperiodic_pincat.mat')

for t = 1:50
    imshow(abs(new(:,:,t)),[]); pause(0.1);
end


% load('invivo_perfusion.mat')
% for t = 1:70
%     imshow(abs(x(:,:,t)),[]); pause(0.1);
% end
% new = x;

%% x-t domain
xidx = 60;
xt = squeeze(new(xidx,:,:));
mshow(xt)

%% x-f domain
xf = ifftshift(fft(fftshift(xt,2),[],2),2);
mshow(xf);

% xtb = ifftshift(ifft(fftshift(xf,2),[],2),2);
% mshow(xtb);
% 
% mm = xt-xtb;
% mshow(mm);
%%
[u,s,v] = svd(xt);

mshow(u*s);