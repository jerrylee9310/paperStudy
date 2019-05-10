function res = directsense(kdata,mask,smap)

[nx,ny,nc] = size(kdata);
res = zeros(nx,ny);
b = sum(conj(smap).*ifft2(mask.*kdata),3);

r = b;
p = r;

% parameters setting
tol = 1e-30;
maxiter = 100;
err0 = r(:)'*r(:);

%% main loop
for i = 1:maxiter
    if rem(i,10) == 0
        disp(['Itr = ' int2str(i)]);
    end
    
    aNu = r(:)'*r(:);
        Sp = smap.*p;
        FSp = mask.*fft2(Sp);
        FhFSp = ifft2(FSp);
        ShFhFSp = sum(conj(smap).*FhFSp,3);
    Ap = ShFhFSp;
    aDe = p(:)'*Ap(:); % p'*A*p
    alpha = aNu / aDe;
    
    res = res + alpha*p;
    r_nxt = r - alpha*Ap; % r - alpha*A*p
    
    beta = (r_nxt(:)'*r_nxt(:)) / (r(:)'*r(:));
    
    p = r_nxt + beta*p;
    
    err = r_nxt(:)'*r_nxt(:);
    if err/err0 < tol
        break
    end
    
    r = r_nxt;
end