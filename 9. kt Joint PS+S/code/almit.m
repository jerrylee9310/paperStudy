function [recon, cst] = almit(ud,params,opt)

    if (nargin == 2)
        opt = 'default';
    end

    % parameter setting
    sz = params.size; np = sz(1); nf = sz(2); nt = sz(3);
    aph = params.alpha;
    L = params.psL;
    
    % find Vt
    [Us, Vt] = getVt(ud,L);
    Vf = 1/sqrt(nt)*fft(Vt,[],2);
       
    cls = @(x,r) reshape(x,np*nf,r);
    clst = @(x,r) reshape(x,np,nf,r);
    
    %% main loop
    
    % init
    G = zeros(np*nf,nt);
    cst = costf(ud,Us,Vt,Vf,G,params);
    bb = cls(ud,nt)*Vt';
    
    % iteration
    fprintf('Main  |  Inner  |  Relative Error\n');
    for i = 1:params.al.maxit
        for j = 1:params.al.maxitin
            % 1st subproblem
            G = shrinkage(Us*Vf,aph);
            % 2nd subproblem
            Us_nxt = cgit(Us,Vt,Vf,G,bb,params,opt);
            
            % stopping criterion
            relchange = norm(Us_nxt(:) - Us(:))/norm(Us(:));
            fprintf(' [%d]  |    %d    |     %f\n',i,j,relchange);
            if relchange < params.al.tol
                break
            end
            Us = Us_nxt;
        end
        
        aph = params.decfac*aph;
        cst = [cst costf(ud,Us,Vt,Vf,G,params)];
        fprintf('\n',i);
    end
    recon = clst(Us*Vt,nt);
end
    

function res = costf(ud,Us,Vt,Vf,G,params)

    [npe,nfe,nfr] = size(ud);

    % constants 
    lmbd = params.sparsity; % regularization weight
    alpha = params.alpha; % Huber function parameter
    mask = params.mask ; % mask
    
    % 1. data fidelity term
    tmp1 = ud - mask.*fft2(reshape(Us*Vt,npe,nfe,nfr));
    fidterm = norm(tmp1(:)).^2;
    
    % 2. regularization term 
    tmp2 = Us*Vf - G;
    regterm =  norm(tmp2,'fro').^2;
    
    % 3. aux term
    auxterm = norm(G(:),1);
    
    res = fidterm + (lmbd/(2*alpha))*regterm + lmbd*auxterm;
end

