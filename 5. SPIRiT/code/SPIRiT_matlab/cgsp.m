function [res, rmse] = cgsp(y,GOP, maxit, tol, x0, option)
%
%  Input:
%		y -	Undersampled k-space data
%		GOP -	the SPIRiT kernel operator obtained by calibration
%		maxit -	Maximum number of iterations
%		tol-	tolerance
%       x0    - Initil value
%
% Outputs:
%		res - Full k-space matrix
%       rmse - root mean square error
%
% original object function : Ax_ = b where A:(G-I)Dc' and b:-(G-I)D'y
% => A'Ax_ = A'b ([A]x_ = bb) 

idx_acq = find(abs(y)>0);
idx_nacq = find(abs(y)==0);

yy = GOP*y; 

if strcmp(option,'lsqr')
    [tmpres,rmse] = cgiteration(-yy,tol,maxit,x0,GOP,idx_acq,idx_nacq); % x0(idx_nacq) = Dcx

    res = y;
    res(idx_nacq) = tmpres;

% Conjugate gradient method
elseif strcmp(option, 'cg')
    [sx,sy,nCoils] = size(x0);
    b = getA(-yy,idx_acq,idx_nacq,sx,sy,nCoils,GOP,'transp');
    
    tmp_Ax = getA(x0(idx_nacq),idx_acq,idx_nacq,sx,sy,nCoils,GOP,'notransp');
    Ax = getA(tmp_Ax,idx_acq,idx_nacq,sx,sy,nCoils,GOP,'transp');

    x = x0;
    r = b-Ax;
    r0_norm = norm(r(:));
    rmse = zeros(maxit+1,1);
    rmse(1) = 1;
    p = r;

    for ii=1:maxit
        fprintf('%d-th iteration\n',ii);
        tmp_Ap = getA(p,idx_acq,idx_nacq,sx,sy,nCoils,GOP,'notransp');
        Ap = getA(tmp_Ap,idx_acq,idx_nacq,sx,sy,nCoils,GOP,'transp');
%         Ap = Ap * r0_norm/norm(Ap);
%         alpha = abs(norm(r(:))/(p(:)'*Ap(:)));
        alpha = (r(:)'*r(:))/(p(:)'*Ap(:));  

        x(idx_nacq) = x(idx_nacq) + alpha*p;
        r_nxt = r - alpha*Ap;

        rmse(ii+1) = norm(r_nxt)/r0_norm;
        if ( norm(r_nxt(:))/r0_norm < tol ) || (ii==maxit)
            tmpres = x;
            break
        end

        beta = (r_nxt(:)'*r_nxt(:))/(r(:)'*r(:));
        p = r_nxt + beta*p;
        r = r_nxt;

    end

    res = y;
    res(idx_nacq) = tmpres(idx_nacq);
    end
    


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