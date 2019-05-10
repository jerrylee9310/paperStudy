function[RR] = compute_RR(params)
%% Load params
nx = params.nx;
ny = params.ny;
nlev = params.Wavelet.nlev;

%% Get the frequency responses of wavelet filters at various levels
[WRA] = compute_SWTFiltResp(params);
DCR = ones(nx,ny) - WRA(:,:,3*nlev+1);

%% Get the frequency response of FD filters
Rh = fft2([-1 1], nx, ny);
Rv = fft2([-1 1]', nx, ny);
Rabs = abs(Rh).^2 + abs(Rv).^2;
% Rhh = fft2([1 -1], nx, ny);
% Rvv = fft2([1 -1]', nx, ny);
% Rabs = Rh.*Rhh + Rv.*Rvv;

%% Compute W^T*W and R^T*R
% RR = Rabs + DCR;
RR = Rabs + 1;
