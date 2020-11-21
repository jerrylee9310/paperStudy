clear all; clc; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [CONTENTS]
%   This code grids a cartesian k-space data into a radial data.
%   It also shows difference image between reference image and
%   grid-processed image.
%   %--- Input ----%
%       Cartesian trajectory k-Space data // [nx, ny, # of coils]
%   %--- Output ---%
%       Reference image (from Cartesian data)
%       Grid-processed image (cart-radial-cart image)
%       Difference image (it reveals artifacts)
%
% [NOTE]
%   - Not sure whether I should use oversampling size deappodisation or
%   just original size deappodisation window.
% 
% [MALFUNCTION]
%   - Severe artifcats on phantom image
%
%                                                  @ Jerry, Dec.28.2018 @
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Brain2D.mat % Cartesian k-Space data
% data_ck = DATA; 

data_ck = fft2c(phantom(256));

% psf = zeros(256);
% psf(128,128) = 1;
% data_ck = fft2c(psf);
%% parameters
fov = size(data_ck,1);
spk = 128; % number of spokes
nc = size(data_ck,3); % number of coils
osf = 2;
kern_wid = 3;
kern_const = 7;
eps = 0.1; % pitty size adjustment for deappodisation

%% setup
zpad = zeros(osf*fov,osf*fov,nc);
data_rk = zeros(fov,spk,nc);
data_rk(fov/2,:,:) = 1;
traj = getTraj(data_rk);
im_size = fov*osf;
ramlak = abs(linspace(-1,1,fov)).';
deap = deappo(kern_const,kern_wid,im_size); % im_size or fov
% deap = deappo(kern_const,kern_wid,fov) + eps; % im_size or fov
% deap = 1;


%% Gridding from Cartesian to Radial Trajectory
% oversampling
if size(deap,1) == fov
    zpad(((im_size-fov)/2 +1):((im_size+fov)/2),((im_size-fov)/2 +1):((im_size+fov)/2),:) = fft2c(data_ck)./deap;
else
    zpad(((im_size-fov)/2 +1):((im_size+fov)/2),((im_size-fov)/2 +1):((im_size+fov)/2),:) = fft2c(data_ck);  
    zpad = zpad./deap;
end
os_ck = ifft2c(zpad);

% Gridding(Cart -> Rad)
for c = 1:nc
    data_rk(:,:,c) = invgridkb(os_ck(:,:,c),traj(:,:,c),fov,osf,kern_wid,kern_const);
end

%% Regridding from Radial to Cartesian Trajectory
% Gridding(Rad -> Cart)
re_ck = zeros(im_size,im_size,nc);
denc = ramlak.*data_rk; % density correction
for c_idx = 1:nc
    re_ck(:,:,c_idx) = gridkb(denc(:,:,c_idx),traj(:,:,c_idx),1,fov,osf,kern_wid,kern_const);
end

re_img = fft2c(re_ck);
if size(deap,1) == fov
    trim_img = re_img(((im_size-fov)/2 +1):((im_size+fov)/2),((im_size-fov)/2 +1):((im_size+fov)/2),:)./deap;
else
    re_img = re_img./deap;
    trim_img = re_img(((im_size-fov)/2 +1):((im_size+fov)/2),((im_size-fov)/2 +1):((im_size+fov)/2),:);
end

%% termination

ss_img = sqrt(sum(abs(trim_img).^2,3));
ref_img = sqrt(sum(abs(fft2c(data_ck)).^2,3));
diff_img = ref_img./max(max(abs(ref_img))) - ss_img./max(max(abs(ss_img)));

% figure();
% subplot(1,3,1); imshow(abs(ref_img./max(max(abs(ref_img)))),[]); title('REFERENCE'); colorbar;
% subplot(1,3,2); imshow(abs(ss_img./max(max(abs(ss_img)))),[]); title('RECONSTRUCTION'); colorbar;
% subplot(1,3,3); imshow(abs(diff_img),[]); title('DIFFERENCE'); colorbar;

% figure(); imshow(abs(ss_img),[0 2e-10]); colorbar; % PSF test
figure(); plot(abs(ss_img(fov/2,:)));
