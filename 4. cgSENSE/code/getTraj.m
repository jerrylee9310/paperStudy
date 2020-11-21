function traj = getTraj(kdata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Contents]
%   It calculate the cartesian point of each radial data point
%   Range of cartesian is from -0.5 to 0.5
%   %--- Input ----%
%       kdata : radial data // [fov, spoke, # of coils]
%   %--- Output ---%
%       traj : real = kx, image = ky
%
%                                                  @ Jerry, Dec.26.2018 @
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fov,spk,nc] = size(kdata);

traj = zeros(size(kdata));

dth = 360/spk; % interval between each spoke
tht = linspace(0,360-dth,spk);

for c_idx =1:nc
    [M,I] = max(kdata(:,:,c_idx));
    peak = mean(I); % calculate the DC point
    radi = ([(1:fov)]-peak)/fov; % radius, -0.5 ~ 0.5
    [x,y]=pol2cart(deg2rad(tht),radi.');
    traj(:,:,c_idx) = x+1j*y;
end

