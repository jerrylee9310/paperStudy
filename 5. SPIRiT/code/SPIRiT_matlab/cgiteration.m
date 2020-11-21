function [x,rmse] = cgiteration(b,tol,maxit,x0,varargin)
%
% Solve "|b-Ax|^2 => A'Ax = A'b (|Ax =|b)" base on 'LSQR' method
%   https://web.stanford.edu/class/cme324/paige-saunders2.pdf
%
% SPIRiT case :
%       A = (G-I)Dc'
%       x = x_tilde(k-space data on unacquired position)
%       b = -(G-I)D'y // y:undersampling data
%
%           => solve A'Ax = A'b // [Dc(G'-I)(G-I)Dc']*[x_tilde] = [Dc(G'-I)(G-I)D'y]
%
%   INPUT:
%       b       : (G-I)D'*y
%       tol     : tolerance
%       maxit   : maximum iteration number
%       x0      : inital guass (zeros)
%       varargin: kernel operator(1), acquired data position(2), non-acquired data position(3)
%

% variable setting
GOP = varargin{1};
[sx,sy,nCoils] = size(x0);
idx_acq = varargin{2};
idx_nacq = varargin{3};
rmse = zeros(maxit+1,1);     % Preallocate vector for norm of residuals

x = x0(idx_nacq); % x_tilde
d = zeros(size(x,1),1);

n2b = norm(b(:));                     % Norm of rhs vector, b
tolb = n2b;                  % Relative tolerance

%% initilisation
u = b - getA(x,idx_acq,idx_nacq,sx,sy,nCoils,GOP,'notransp'); % u = b-Ax // residual
beta = norm(u(:));
u = u / beta;
normr = beta;

% iteration parameters
rmse(1) = normr;             % resvec(1,1) = norm(b-A*x0)
phibar = beta;
c = 1;
s = 0;

v = getA(u,idx_acq,idx_nacq,sx,sy,nCoils,GOP,'transp'); % v = A'b-A'Ax = |b - |Ax // LS residual
alpha = norm(v);
v = v / alpha;

%% main loop
for ii = 1 : maxit
    fprintf('%d-th iteration\n',ii);
    
    z = v;
    
    % beta(t+1)*u(t+1) = Av(t) - alpha(t)*u(t)
    u = getA(z,idx_acq,idx_nacq,sx,sy,nCoils,GOP,'notransp') - alpha * u;
    beta = norm(u(:));
    u = u / beta;
    
    % parameter updates
    thet = - s * alpha;
    rhot = c * alpha;
    rho = sqrt(rhot^2 + beta^2);
    c = rhot / rho;
    s = - beta / rho;
    phi = c * phibar;
    phibar = s * phibar;
    d = (z - thet * d) / rho;
    
    x = x + phi * d;
    normr = abs(s) * normr;
    rmse(ii+1) = normr;
    
    % stop condition
    if normr < tol*tolb
        break
    end
    
    vt = getA(u,idx_acq,idx_nacq,sx,sy,nCoils,GOP,'transp');

    v = vt - beta * v;
    alpha = norm(v);
    v = v / alpha;
    
end                      

rmse = rmse/tolb;



function res = getA(x,idx_acq,idx_nacq,sx,sy,nCoils,GOP,flag)

if strcmp(flag,'notransp')
    tmp = zeros(sx,sy,nCoils);
    tmp(idx_nacq) = x;
    res = GOP*tmp;
elseif strcmp(flag,'transp')
    tmp = reshape(x,sx,sy,nCoils);
    res = GOP'*tmp;
    res = res(idx_nacq);
end



