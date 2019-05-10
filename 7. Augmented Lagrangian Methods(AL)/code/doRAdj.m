function [Rz] = doRAdj(z, params)
% Do R*z, where R is a linear operator that signifies a wavelet transform or finite differences

%% Load parameters
nlev = params.Wavelet.nlev;
lor = params.Wavelet.lor;
hir = params.Wavelet.hir;

Operator = params.Operator;

% %% Compute v (Refer notes)       
% [h, v] = trfd(z(:, :, 1), z(:, :, 2));
% Rz = h + v;
% z1 = z(:, :, 3:3*nlev+3);
% Rz = Rz + myiswt2(z1, lor, hir);
%             

%% Compute v (Refer notes)
switch(Operator)
    case{'FD'} % Finite difference
        [h, v] = trfd(z(:, :, 1), z(:, :, 2));
        Rz = h + v;
        
    case{'W'} % Wavelets
        z1 = z;
        Rz = myiswt2(z1, lor, hir);
        
    case{'WFD'}
        [h, v] = trfd(z(:, :, 1), z(:, :, 2));
        Rz = h + v;
        z1 = z(:, :, 3:3*nlev+3);
        Rz = Rz + myiswt2(z1, lor, hir);
end