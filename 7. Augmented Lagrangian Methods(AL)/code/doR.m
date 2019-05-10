function [Rz] = doR(z, params)
% Do R*z, where R is a linear operator that signifies a wavelet transform or finite differences

%% Load parameters
nlev = params.Wavelet.nlev;
lod = params.Wavelet.lod;
hid = params.Wavelet.hid;

Operator = params.Operator;

% %% Compute v (Refer notes)
% 
%         % And finite difference
% [h, v] = fd(z);
% Rz(:, :, 1) = h;
% Rz(:, :, 2) = v;
% 
% 
% SWC = myswt2(z, nlev, lod, hid);
% Rz(:, :, 3:3*nlev+3) = SWC;


%% Compute v (Refer notes)
switch(Operator)
    case{'FD'} % Finite difference
        [h, v] = fd(z);
        Rz(:, :, 1) = h;
        Rz(:, :, 2) = v;
        
    case{'W'} % Wavelets
        SWC = myswt2(z, nlev, lod, hid);
        Rz = SWC;
            
    case{'WFD'}
        % And finite difference
        [h, v] = fd(z);
        Rz(:, :, 1) = h;
        Rz(:, :, 2) = v;
        % wavelet
        SWC = myswt2(z, nlev, lod, hid);
        Rz(:, :, 3:3*nlev+3) = SWC;
end