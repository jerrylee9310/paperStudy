function out = ktRPCA_ADMM(y,opts)
% k-t RPCA using alternating direction method of multipliers
%
% Solve
% min_{L,S} 0.5||E(L+S)-y||^2_2 + mu1 ||L||_* + mu2 ||T(S)||_1 
%
% Ref:
% Tremoulheac et al, Dynamic MR image reconstruction--separation from 
% undersampled (k,t)-space via low-rank plus sparse prior, TMI 2014
%
%
% Benjamin Trémoulhéac
% b.tremoulheac@cs.ucl.ac.uk


mu1 = opts.mu1; % nuclear norm
mu2 = opts.mu2; % L1 norm

tol = opts.tol;
maxIter = opts.maxIter;
E = opts.E;
Eh = opts.Eh;
EhE = opts.EhE;
T = opts.T;
Tt = opts.Tt;
nx = opts.datasize(1);
ny = opts.datasize(2);
nt = opts.datasize(3);
quiet = opts.quiet;

delta1 = 1;
delta2 = 1;

X = Eh(y);
L = X;
S = zeros(nx*ny*nt,1);
Z1 = zeros(nx*ny*nt,1);
Z2 = zeros(nx*ny*nt,1); 

it = 0;

fprintf('k-t RPCA reconstruction\n');

if quiet == false
    fprintf('It.\tData fidelity\tNuclear norm\tl1 norm\t\tObjective\tStop criterion\n');
end

converged = false;

while (~converged && it < maxIter)
    Xp = X;
    P = svt(reshape(X(:)-S(:),[nx*ny nt]) + (1/delta1)*reshape(Z1,[nx*ny nt]),mu1/delta1);
    Q = shrink(T(X-L) + (1/delta2)*Z2(:), mu2/delta2);
    
    bl = Eh(y) + delta1*P(:) - Z1(:) - EhE(S);
    [L,~] = cgs(@(z) EhE(z(:)) + delta1*z(:),bl,1e-5,100);
    
    bs = Eh(y) + Tt(delta2*Q(:) - Z2(:)) - EhE(L);
    [S,~] = cgs(@(z) EhE(z(:)) + delta2*z(:),bs,1e-5,100);

    Z1 = Z1+ delta1*(L-P(:));
    TS = T(S);
    Z2 = Z2 + delta2*(TS-Q);

    X = L+S;
    it = it+1;
    
    out.stopCriterion(it) = norm(X(:)-Xp(:))/norm(Xp(:));
    
    ELSy = E(X)-y;
    out.fval1(it) = 0.5*norm(ELSy(:),2)^2;
    out.fval2(it) = mu1*norm_nuc(reshape(L,[nx*ny nt]));
    out.fval3(it) = mu2*norm(TS(:),1);
    out.fval(it) = out.fval1(it) + out.fval2(it)+ out.fval3(it);

    if quiet == false
        fprintf('%g\t%e\t%e\t%e\t%e\t%e\n',it,out.fval1(it),out.fval2(it),out.fval3(it),out.fval(it), out.stopCriterion(it));
    end
    
    if it == maxIter
        fprintf('Stopped because maximum of iterations (%g) reached\n',maxIter);
    end
    
    if out.stopCriterion(it) <= tol
        converged = true;
        if quiet == false
            fprintf('Stopped because TOL=%.1e achieved\n',tol);
        end
    end
    
end
out.recon = X;
out.L = L;
out.S = S;
out.nb_iter = it;
end

function B = svt(A,gamma)
[U,S,V] = svd(A,0);
Z = diag(shrink(diag(S),gamma));
B = U*Z*V';
end

function A = shrink(B,gamma)
A = sign(B).*max(abs(B)-gamma,0);
end

function B = norm_nuc(A)
B = sum(svd(A,'econ'));
end
