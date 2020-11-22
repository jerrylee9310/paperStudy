function res = cgit(Us,Vt,Vf,G,bb,params,opt)

    sz = params.size; np = sz(1); nf = sz(2); nt = sz(3);
    mu = params.sparsity;
    aph = params.alpha;
    L = params.psL;
    mask = params.mask;
    tol = 1e-10;
    iter = 200;
    
    cls = @(x,t) reshape(x,np*nf,t);
    clst = @(x,t) reshape(x,np,nf,t);
   
    if opt == 'Cart'
        % init
        Uk = cls(fft2(clst(Us,L)),L);
        rhs = bb + (mu/(2*aph))*cls(fft2(clst(G*Vf',L)),L);

        % iteration
        [Uk_nxt,~] = pcg(@(x) (normal_oper(x,Vt,mask) + (mu/(2*aph))*x), rhs(:), tol, iter, [], [], Uk(:));
        Uk_nxt = cls(Uk_nxt,L);
        res = cls(ifft2(clst(Uk_nxt,L)),L);
    else
        % init
        rhs = bb + (mu/(2*aph))*cls(fft2(clst(G*Vf',L)),L);

        % iteration
        [Uk_nxt,~] = pcg(@(x) (normal_operc(x,Vt,mask) + (mu/(2*aph))*x), rhs(:), tol, iter, [], [], Us(:));
        Uk_nxt = cls(Uk_nxt,L);
        res = cls(ifft2(clst(Uk_nxt,L)),L);
    end
    
end