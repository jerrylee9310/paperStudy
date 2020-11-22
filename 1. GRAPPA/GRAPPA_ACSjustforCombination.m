
clear all; close all; clc;

%% load rawdata
load Brain2D.mat
% CoilData = fftshift(fft2(fftshift(CoilIm)));
raw = DATA;
%% parameter
Acc_factor = 2;
Acs_num = 35; % How to modify the code when it has more than one ACS line?
Nb = 4;
m = 3;

% Nb : the number of blocks using in calibration & reconstruction



%% get size of each dimension
[Ny, Nx, Nc] = size(raw);
% Ny : y-dimension
% Nx : x-dimension
% Nc : coil-dimension
%
% NOTE : In matlab the first index stands for row    (which is vertical)
%                    second index stands for column (which is horizontal)

%% Making each coil image by inverse Fourier transform
coil_img = ifftshift(ifft2(ifftshift(raw)));
% NOTE : ifft2 and ifftshift does not affect coil dimension
% for more information, type "doc ifft2" and "doc ifftshift" on Command
% Window

subplot(1,2,1);
imshow(abs(sum(coil_img(:,:,:),3)),[])

%% Displaying each coil image
% figure('Name', 'Each coil image'),
% for c = 1:Nc
%     subplot(2,ceil(Nc/2),c),
%     
%     imshow(abs(coil_img(:,:,c)),[])
% end

%% Undersampling the data
under_raw = zeros(Ny,Nx,Nc);
under_array = 1:Acc_factor:Ny;
len = length(under_array);
Acs_index = len/2; % ACS line is acquired at the center of raw data
Acs = under_array(Acs_index)+1; % NB in the case you want to shift the ACS, you can change the integar smaller than abs(Acc_factor)

under_raw(under_array,:,:) = raw(under_array,:,:);
% Acs lines
upper = Acs:Acc_factor:Acs+floor(Acs_num/2)*Acc_factor;
lower = Acs:-Acc_factor:Acs-round(Acs_num/2)*Acc_factor+1;
under_raw(upper,:,:) = raw(upper,:,:);
under_raw(lower,:,:) = raw(lower,:,:);

%% Calibration

% initailise
A = zeros(Nx,Nc*Nb);
B = zeros(Nx,Nc);
k_y = Acs + m; % sliding block approach (starting point)

% Acquired line data
for l = 1:Nc
    for b = 0:(Nb-1)
        A(:,Nb*(l-1)+b+1) = under_raw(k_y-b*Acc_factor,:,l)'; % sliding block approach (overall shape)
    end
end

% Acs line data
for l = 1:Nc
    B(:,l) = under_raw(Acs,:,l)';
end

% Get weighting
weighting = pinv(A)*B;
    
%% Reconstruction
n = (Nb-1)*Acc_factor-m-1;
rec_Ny = Ny+n+m;
recon_raw = zeros(rec_Ny,Nx,Nc);

% High frequency handling
recon_raw(n+1:Ny+n,:,:) = under_raw(:,:,:);
recon_raw(1:n,:,:) = under_raw(Ny-(n-1):Ny,:,:);
recon_raw(Ny+n+1:rec_Ny,:,:) = under_raw(1:m,:,:);
%%
for i = 1:2:lower(length(lower))-3
    o = permute(recon_raw(i:Acc_factor:i+Acc_factor*(Nb-1),:,:),[2 1 3]);
    C = o(:,:,1);
    for j=2:Nc
        C = cat(2,C,o(:,:,j));
    end   
    recon_raw(i+n+1,:,:) = C*weighting;
end

for i = upper(length(upper))+1:2:Ny
    o = permute(recon_raw(i:Acc_factor:i+Acc_factor*(Nb-1),:,:),[2 1 3]);
    C = o(:,:,1);
    for j=2:Nc
        C = cat(2,C,o(:,:,j));
    end   
    recon_raw(i+n+1,:,:) = C*weighting;
end

recon_data = recon_raw(n+1:n+Ny,:,:);
    
%% Display reconstruction image
recon_img = ifftshift(ifft2(ifftshift(recon_data)));

squared_img = abs(recon_img).^2;
sum_of_squared_img = sum(squared_img, 3);
sos_img = sqrt(sum_of_squared_img);

subplot(1,2,2);
imshow(sos_img,[]);