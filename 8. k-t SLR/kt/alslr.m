function [res,cost] = alslr(Fspat,Spar,udata,params)

    % regularization weight
    mu1 = params.reg.nuclear;
    mu2 = params.reg.sparsity;

    % penalty weight
    beta1 = params.constr.casorati;
    beta2 = params.constr.sparsity;

    inc1 = params.al.incfact.casorati;
    inc2 = params.al.incfact.sparsity;


    %% init
    res = zeros(size(udata)); % reconstructed matrix: T
    R = res; % 1st splitted variable: T = R
    S = Spar*res; % 2nd splitted variable: S = Q(T)

    % Lagrage multiplier
    eta1 = zeros(size(R)); % for 'R-T' constraint 
    eta2 = zeros(size(S)); % for 'S-Q(T)' constraint

    %% AL
    cost = costf(Fspat,Spar,udata,res,R,S,eta1,eta2,params);

    for i = 1:params.al.maxit
        fprintf('%d\n',i);
        for j = 1:params.al.maxitin
            % 1st subproblem
            S = shrinkage(Spar*res + eta2/beta2,mu2/beta2,'mul');
            
            % 2nd subproblem
            R = shrinkage(res+eta1/beta1,mu1/beta1,'sig');

            % 3rd subproblem
            res = cgit(Fspat,Spar,res,udata,R,S,eta1,eta2,params);

%             eta1 = eta1 - beta1*(R - res);
%             eta2 = eta2 - beta2*(S - Spar*res);
        end

        beta1 = inc1*beta1;
        beta2 = inc2*beta2;

        % break condition
        cost = [cost costf(Fspat,Spar,udata,res,R,S,eta1,eta2,params)];
        if cost(end)/cost(1) < params.al.tol
            break
        end
        
%         imshow(abs(res(:,:,1)),[]);
    end
end


function res = costf(Fspat,Spar,udata,res,R,S,eta1,eta2,params)

    % regularization weight
    lmbd1 = params.reg.nuclear;
    lmbd2 = params.reg.sparsity;

    % penalty weight
    beta1 = params.constr.casorati;
    beta2 = params.constr.sparsity;
    
    % 1. data fidelity term
    tmp1 = Fspat*res - udata;
    fidterm = tmp1(:)'*tmp1(:);
    
    % 2. regularization term
    tmp2 = R - (res - eta1/beta1);
    regterm = nnorm(R) + 0.5*beta1*(tmp2(:)'*tmp2(:));
    
    % 3. prior term
    tmp3 = S - (Spar*res + eta2/beta2);
    prterm = getPsi(S,params) + 0.5*beta2*(tmp3(:)'*tmp3(:));
    
    res = fidterm + lmbd1*regterm + lmbd2*prterm;
end



    
    
        