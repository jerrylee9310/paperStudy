function res = freqrespver(Fspat,Spar,res,udata,R,S,eta1,eta2,params)

[ny,nx,nt] = size(res);

mu1 = params.reg.nuclear;
mu2 = params.reg.sparsity;

beta1 = params.constr.casorati;
beta2 = params.constr.sparsity;

% solve Ax=b using lsqr 
b = 2*(Fspat'*udata) + mu1*beta1*R - mu1*eta1 + mu2*beta2*(Spar'*S) - mu2*(Spar'*eta2);

% Inverse of some matrices required for solving sub-problems
Rh = fft2([-1 1], ny, nx);
Rv = fft2([-1 1]', ny, nx);
Rabs = abs(Rh).^2 + abs(Rv).^2;


Dt = imfilter(x,reshape([-1 1],1,1,2),'circular');
iDt = cat(3,Dt(:,:,end),Dt);
iDt = imfilter(iDt,reshape([1 -1],1,1,2), 'circular');
iDt = iDt(:,:,1:end-1);
            
den = 2*params.mask + mu1*beta1 + mu2*beta2*Rabs;

tmp = fft2(b)./den;
res = ifft2(tmp);