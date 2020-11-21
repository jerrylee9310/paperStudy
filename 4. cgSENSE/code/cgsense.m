clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [REFERENCE]
%   Advances in Snsitivity Encoding With Arbirtary k-Space Trajectories,
%   Pruessmann, MRM, 2001
%
% [CONTENTS]
%   This code reconstruct a radial data using conjugate gradient method
%   %--- Input ----%
%       radial k-space data // [fov, spoke, # of coils]
%   %--- Output ---%
%       reconstruction image
%
% [NOTE]
%   - Keyword : "Congugate Gradient", "Radial k-space data", "SENSE"
%   - I only test this code on radial data, so I'm not sure bout other
%     trajectoreis LMAO.
%   - If you want to generate radial data from cartesian data, use
%     'Cart2Rad.m'.
%   - No intensity correction included (disadvantage on computational efficiency)
% 
% [MALFUNCTION]
%   - deappodization
%   - only work on brain2RK // numerical phantom, brainRK -> X
%
%                                                  @ Jerry, Dec.26.2018 @
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load brain2RK.mat;
kdata = brain2RK;

[fov,spk,nc] = size(kdata);
%% manipulable variables

% BASIC VARS
acc = 4; % reduction factor

% GRIDDING VARS
osf = 2; % oversampling factor
kern_wid = 3; % kernel width
kern_const = 7; % kernel constant for kaisser-bessel window
eps = 0.1; % pitty size adjustment for deappodisation

% CG LOOP VARS
stp = 1e-4; % stop condition(accuracy)
max_iter = 10; % maximum iteration number

%% parameters setting

% relative parameters
im_size = osf*fov;
traj = getTraj(kdata); % radial data's trajectory, -0.5 ~ 0.5
ramlak = abs(linspace(-1,1,fov)).'; % Ram-Lak filter for density correction
lowfrq_wid = 24; % low frequency width for coil sensitivity
zpad = zeros(im_size,im_size); % for crop and 2X oversampling
zpad(1+((osf-1)*fov)/2:((osf+1)*fov)/2, 1+((osf-1)*fov)/2:((osf+1)*fov)/2) = 1;
% deap = 1; 
deap = zeros(im_size);
deap(1+((osf-1)*fov)/2:((osf+1)*fov)/2, 1+((osf-1)*fov)/2:((osf+1)*fov)/2) = deappo(kern_const,kern_wid,fov);
deap = deap+eps;

% coil sensitivity
coil_sen = zeros(im_size,im_size,nc);
cs_mask = zeros(fov,spk);
cs_mask((fov-lowfrq_wid)/2:(fov+lowfrq_wid)/2 -1,:) = repmat(hanning(lowfrq_wid),[1,spk]);
for c_idx = 1:nc
    cs_tmp = ramlak.*(cs_mask.*kdata(:,:,c_idx));
    coil_sen(:,:,c_idx) = fft2c(gridkb(cs_tmp,traj(:,:,c_idx),1,fov,osf,kern_wid,kern_const));
end
coil_sen = coil_sen./sqrt(sum(abs(coil_sen).^2,3));

% reference image
load brain2RF;
tmp = zeros(im_size,im_size);
tmp(1+((osf-1)*fov)/2:((osf+1)*fov)/2, 1+((osf-1)*fov)/2:((osf+1)*fov)/2) = brain2RF;
ref_img = tmp;
% coil_img = zeros(im_size,im_size,nc); 
% for c_idx = 1:nc
%     cs_tmp = ramlak.*kdata(:,:,c_idx);
%     coil_img(:,:,c_idx) = fft2c(gridkb(cs_tmp,traj(:,:,c_idx),1,fov,osf,kern_wid,kern_const));
% end
% ref_img = zpad.*sqrt(sum(abs(coil_img).^2,3));

%% Initialisation

% undersampling
red_k = kdata(:,1:acc:end,:);

% start image
strt_cim = zeros(im_size,im_size,nc);
for c_idx = 1:nc
    % density correction
    im_tmp = ramlak.*red_k(:,:,c_idx); 
    % grdding
    strt_cim(:,:,c_idx) = (fft2c(gridkb(im_tmp,traj(:,1:acc:end,c_idx),1,fov,osf,kern_wid,kern_const)))./deap; 
    % multiply conjugate coilsensitivity
    strt_cim(:,:,c_idx) = conj(coil_sen(:,:,c_idx)).*strt_cim(:,:,c_idx);
end
strt_im = sum(strt_cim,3);
strt_im = zpad.*strt_im;

% storage for reconstructed image at each loop
cg_recon = zeros(im_size,im_size,max_iter); 
dif_img = zeros(im_size,im_size,max_iter); 
rmse_lv = zeros(1,max_iter);
%% main(cg loop)

% Loop init
srch_dir = strt_im(:); % search direction
res = strt_im(:); % residuum
recon = zeros(size(res)); % reconstructed image

% CG loop
for it = 1:max_iter
    
    % multiply coisensitivity
    FT2img = coil_sen.*reshape(srch_dir,im_size,im_size)./abs(deap);
    
    % FT2 // regridding (inverse gridding)
    FT2ck = ifft2c(FT2img);
    FT2rk = zeros(fov,spk,nc);
    for c = 1:nc
        FT2rk(:,:,c) = invgridkb(FT2ck(:,:,c),traj(:,:,c),fov,osf,kern_wid,kern_const);
    end
    
    % masking
    maskrk = FT2rk(:,1:acc:end,:);
    
    % FT1 // gridding (forward gridding)
    FT1 = zeros(im_size,im_size,nc);
    for c=1:nc
        im_tmp = ramlak.*maskrk(:,:,c); % density correction
        FT1(:,:,c) = gridkb(im_tmp,traj(:,1:acc:end,c),1,fov,osf,kern_wid,kern_const)./deap;
    end
    
    % multiply conjugate coil sensitivity
    cgim = conj(coil_sen).*fft2c(FT1);
    
    % summation and trimming
    sum_cgim = zpad.*sum(cgim,3);
    q = sum_cgim(:);
    
    % cg parameters update
    alpha = (res'*res)/(srch_dir'*q);
    recon = recon + alpha*srch_dir;
    r_nxt = res - alpha*q;
        
    beta = (r_nxt'*r_nxt)/(res'*res);
    srch_dir = r_nxt + beta*srch_dir;
    res = r_nxt;
    
    dlt = (res'*res)/((strt_im(:))'*(strt_im(:))); % accuracy
    
    % display part
    disp([num2str(it) 'th iteration']) 
    disp(['convergence =' num2str(dlt)])
    rmse_lv(it) = dlt;
    cg_recon(:,:,it) = reshape(recon,im_size,im_size);
    dif_img(:,:,it) = ref_img./max(max(abs(ref_img))) - reshape(recon,im_size,im_size)./max(abs(recon));
    
    subplot(1,2,1); imshow(abs(cg_recon(:,:,it)),[]); drawnow;
    subplot(1,2,2); imshow(abs(dif_img(:,:,it)),[]); drawnow;
    
    
    cg_recon(:,:,it) = reshape(recon,im_size,im_size);
    dif_img(:,:,it) = ref_img./max(max(abs(ref_img))) - reshape(recon,im_size,im_size)./max(abs(recon));
    
    % stop criterion
    if dlt < stp
        break;
        
    end   
end
%% temination
mshow(cg_recon(1+((osf-1)*fov)/2:((osf+1)*fov)/2, 1+((osf-1)*fov)/2:((osf+1)*fov)/2,:),2,5,10);
mshow(dif_img(1+((osf-1)*fov)/2:((osf+1)*fov)/2, 1+((osf-1)*fov)/2:((osf+1)*fov)/2,:),2,5,10);