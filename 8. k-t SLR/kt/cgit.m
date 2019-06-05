function res = cgit(Fspat,Spar,x,udata,R,S,eta1,eta2,params)
[ny,nx,nt] = size(x);

beta1 = params.constr.casorati;
beta2 = params.constr.sparsity;

tol = 1e-20;
maxit = 30;

% solve Ax=b using lsqr 
b = 2*(Fspat'*udata) + (beta1*R - eta1) + Spar'*(beta2*S - eta2);
res = lsqr(@(x,tflag)afun(x,tflag),b(:),tol,maxit,[],[],x(:));
    
    % A = function below
    function res = afun(input,tflag)       
        tmp = reshape(input,ny,nx,nt);
        res = 2*(Fspat'*(Fspat*tmp)) + beta1*tmp + beta2*(Spar'*(Spar*tmp));
        res = res(:);
    end

res = reshape(res,ny,nx,nt);
end



%% CG iteration (line by line code)

% res = x;
% 
% rhs = 2*(Fspat'*udata) + mu1*beta1*R - mu1*eta1 + mu2*beta2*(Spar'*S) - mu2*(Spar'*eta2);
% Ax = 2*(Fspat'*(Fspat*res)) + mu1*beta1*res + mu2*beta2*(Spar'*(Spar*res));
% 
% r = rhs - Ax;
% p = r;
% err0 = r(:)'*r(:);
% 
% for i = 1:maxit
%     fprintf('%d',i');
%     aNu = r(:)'*r(:);
%     Ap = 2*(Fspat'*(Fspat*p)) + mu1*beta1*p + mu2*beta2*(Spar'*(Spar*p));
%     aDe = abs(p(:)'*Ap(:)); % p'*A*p
%     alpha = real(aNu) / real(aDe);
% 
%     res = res + alpha*p;
%     r_nxt = r - alpha*Ap; % r - alpha*A*p
% 
%     beta = (r_nxt(:)'*r_nxt(:)) / (r(:)'*r(:));
% 
%     p = r_nxt + beta*p;
% 
%     err = r_nxt(:)'*r_nxt(:);
%     if err/err0 < tol
%         break
%     end
% 
%     r = r_nxt;
% end
% fprintf('\n');
% end
