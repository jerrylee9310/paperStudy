function res = ktRPCA(d,params)

    nx = params.dsize(1); ny = params.dsize(2); nt = params.dsize(3);

    mu = params.mu;
    lmbd = params.lambda;
    del1 = params.penalty_LR;
    del2 = params.penalty_SP;

    tol = params.tol;
    maxit = params.maxit;
    cgtol = params.cgtol;
    cgit = params.cgit;

    E = params.E;
    Eh = params.Eh;
    EhE = params.EhE;
    Ft = params.Ft;
    Fth = params.Fth;


    %% iteration

    % initialisation
    X = Eh(d);
    L = X; % X = L+S
    S = zeros(nx*ny*nt,1); % X = L+S, Q = Ft(S);
    Z1 = zeros(nx*ny*nt,1); % largrangian multiplier for 'L = P' constraint
    Z2 = zeros(nx*ny*nt,1); % largrangian multiplier for 'Q = Ft(S)' constraint

    % main loop
    it = 0;
    converged = false;

    while (~converged && it < maxit)
        fprintf('%d  ',it); 
        Xp = X;
        % 1st subproblem
        P = shrinkage(reshape(X(:)-S(:),[nx*ny nt]) + reshape(Z1,[nx*ny nt]),mu/del1,'singular');
        
        % 2nd subproblem
        Q = shrinkage(Ft(X-L) + (1/del2)*Z2(:), (mu*lmbd)/del2,'1D');
       
        % 3rd subproblem
        rhs_L = Eh(d) + del1*P(:) - Z1(:) - EhE(S);
        [L,~] = cgs(@(x) EhE(x(:)) + del1*x(:),rhs_L,cgtol,cgit);
        
        % 4th subproblem
        rhs_S = Eh(d) + Fth(del2*Q(:) - Z2(:)) - EhE(L);
        [S,~] = cgs(@(x) EhE(x(:)) + del2*x(:),rhs_S,cgtol,cgit);
        
        % Lagrangain multiplier update
        Z1 = Z1 + del1*(L - P(:));
        Z2 = Z2 + del2*(Ft(S) - Q);

        X = L+S;
        it = it+1;
        stopCondition(it) = norm(X(:)-Xp(:))/norm(Xp(:));
        fprintf('%f\n',stopCondition(it));
        if stopCondition(it) <= tol
            converged = true;
        end
    end
    
    res = reshape(L+S,[nx ny nt]);
end
