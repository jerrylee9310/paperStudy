%% Features %%
%

%%
clear all; close all; clc;

%% load rawdata
load matlab2.mat
CoilData = fftshift(fft2(fftshift(CoilIm)));
raw = CoilData;

%% parameter
Acc_factor = 2;
Acs_num = 30; 
Nb = 4;
m = 2;
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

subplot(2,2,1);
imshow(sqrt(sum(abs(coil_img).^2,3)),[])
title('Reference');

%% Displaying each coil image
% figure('Name', 'Each coil image'),
% for c = 1:Nc
%     subplot(2,ceil(Nc/2),c),
%     
%     imshow(abs(coil_img(:,:,c)),[])
% end

%% Undersampling the data
% initialise
under_raw = zeros(Ny,Nx,Nc);
under_acquired = 1:Acc_factor:Ny;
under_unaquired = setdiff(1:Ny,under_acquired);

% Acs line shape
% 1. alternative(right-left sequence)
% Acs_array = under_unaquired(round(length(under_unaquired)/2)-floor(Acs_num/2):round(length(under_unaquired)/2)+round(Acs_num/2)-1);

% 2. block-priority
num_of_cali_block = fix(Acs_num/(Acc_factor-1));
Acs_start_index = fix(length(under_acquired)/2)-ceil(num_of_cali_block/2);
Acs_array_cali = under_acquired(Acs_start_index)+1:under_acquired(Acs_start_index+1)-1;

for i = 1:num_of_cali_block-1
    Acs_array_cali = [Acs_array_cali,under_acquired(Acs_start_index+i)+1:under_acquired(Acs_start_index+i+1)-1];
end

if rem(Acs_num,(Acc_factor-1)) == 0
    Acs_array_tot = Acs_array_cali;
else
    Acs_array_tot = [Acs_array_cali,Acs_array_cali(end)+1:Acs_array_cali(end)+rem(Acs_num,(Acc_factor-1))+1];
end

% undersampling
under_raw(under_acquired,:,:) = raw(under_acquired,:,:);

% picturing undersampled image
u_coil_img = ifftshift(ifft2(ifftshift(under_raw)));
subplot(2,2,3);
imshow(sqrt(sum(abs(u_coil_img).^2,3)),[])
title('undersampled image without ACS');

% insert ACS lines
under_raw(Acs_array_tot,:,:) = raw(Acs_array_tot,:,:);

% picturing undersampled image with ACSs
under_coil_img = ifftshift(ifft2(ifftshift(under_raw)));
subplot(2,2,4);
imshow(sqrt(sum(abs(under_coil_img).^2,3)),[])
title('undersampled image with ACS');

%% Calibration
% initialise
A = zeros(Nx*num_of_cali_block,Nc*Nb);
B = zeros(Nx*num_of_cali_block,Nc*(Acc_factor-1));

% Acquired line data(A)
for stack = 1:num_of_cali_block
    % block is start from acquired line(acquired line + (Acc_factor-1) unacquired line)
    acs_belongs_block_acq = under_acquired(Acs_start_index)+Acc_factor*(stack-1);
    k_y = acs_belongs_block_acq+m*Acc_factor;
    for l = 1:Nc
        for b = 0:(Nb-1)
            % sliding block approach (overall shape) 
            A((stack-1)*Nx+1:stack*Nx,Nb*(l-1)+b+1) = under_raw(k_y-b*Acc_factor,:,l)'; 
        end
    end
end

% Acs line data(B)
for stack = 1:num_of_cali_block
    for unacq_seq = 1:Acc_factor-1
        Acs = Acs_array_cali((Acc_factor-1)*(stack-1)+unacq_seq);
        for l = 1:Nc
            B((stack-1)*Nx+1:stack*Nx,(unacq_seq-1)*Nc+l) = under_raw(Acs,:,l)';
        end
    end
end

% Get weighting
weighting = pinv(A)*B;
    
%% Reconstruction
% initialise
% vertical appending for high frequency handling
top_block = Nb-m-1;
if top_block<0
    top_block = 0;
end
top_appending = top_block*Acc_factor;

bottom_block = m;
if bottom_block<0
    bottom_block = 0;
end
bottom_appending = bottom_block*Acc_factor;

last_line_trimer = (Acc_factor-1) - (Ny - under_acquired(length(under_acquired)));
rec_Ny = top_appending + Ny + last_line_trimer + bottom_appending;
recon_raw = zeros(rec_Ny,Nx,Nc);

% High frequency handling
recon_raw(top_appending+1:Ny+top_appending,:,:) = under_raw(:,:,:);
recon_raw(1:top_appending,:,:) = recon_raw(top_appending+under_acquired(length(under_acquired)-(top_block-1)):top_appending+under_acquired(length(under_acquired)-(top_block-1))+top_appending-1,:,:);
recon_raw(rec_Ny-bottom_appending+1:rec_Ny,:,:) = recon_raw(top_appending+1:top_appending+bottom_appending,:,:);

for i = 1:Acc_factor:Ny
    o = permute(recon_raw(i+Acc_factor*(Nb-1):-Acc_factor:i,:,:),[2 1 3]); % i:Acc_factor:i+Acc_factor*(Nb-1) gives better result
    A_recon = o(:,:,1);
    for j=2:Nc
        A_recon = cat(2,A_recon,o(:,:,j));
    end
    B_recon = A_recon*weighting;
    for recon_index = 1:Acc_factor-1
        recon_raw(top_appending+i+recon_index,:,:) = B_recon(:,(recon_index-1)*Nc+1:recon_index*Nc);
    end
end

% Using ACS line data before combine the image
recon_raw(Acs_array_tot+top_appending,:,:) = under_raw(Acs_array_tot,:,:); 
recon_data = recon_raw(top_appending+1:top_appending+Ny,:,:);
    
%% Display reconstruction image
recon_img = ifftshift(ifft2(ifftshift(recon_data)));

squared_img = abs(recon_img).^2;
sum_of_squared_img = sum(squared_img, 3);
sos_img = sqrt(sum_of_squared_img);

subplot(2,2,2);
imshow(sos_img,[]);
title('Reconstruced image');