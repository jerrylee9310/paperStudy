function [x] = SENSERecon_ALP1(Fhd, eta, u, x, params)
% u: [u1]  = [Rx]
%%
mask = params.mask;
smap = params.smap;

maxitr = params.maxitr;
maxitr_in = params.maxitr_in;

mu = params.AL.mu;
u1 = u;

itr = 1;

while((itr <= maxitr))
    itr_in = 1; % Inner iterations for minimization over z's
        disp(['ItrALP2 = ' int2str(itr)]);
    
    while(itr_in <= maxitr_in)
        
        % Compute u0 (Refer revised manuscript)
            ShFhd = sum(conj(smap) .* Fhd, 3);
            Rhu1eta = doRAdj(u1-eta,params);
        rhsz = ShFhd + mu*Rhu1eta;
            ShFhFSx = sum(conj(smap).*ifft2(mask.*fft2(smap.*x)),3);
            RhRx = doRAdj(doR(x,params),params);
        Ax = ShFhFSx + mu*RhRx;
        r = rhsz - Ax;
        p = r;
        err0 = r(:)'*r(:);

        %% main loop
        for i = 1:100
            aNu = r(:)'*r(:);
                ShFhFSp = sum(conj(smap).*ifft2(mask.*fft2(smap.*p)),3);
                RhRp = doRAdj(doR(p,params),params);
            Ap = ShFhFSp + mu*RhRp;
            aDe = abs(p(:)'*Ap(:)); % p'*A*p
            alpha = aNu / aDe;

            x = x + alpha*p;
            r_nxt = r - alpha*Ap; % r - alpha*A*p

            beta = (r_nxt(:)'*r_nxt(:)) / (r(:)'*r(:));

            p = r_nxt + beta*p;

            err = r_nxt(:)'*r_nxt(:);
            if err/err0 < 1e-30
                break
            end

            r = r_nxt;
        end
        
        % Compute u1 (Refer revised manuscript)
        RW = doR(x, params); % First do Ru2
        Rw = RW + eta;
        u1 = thresholdRw(Rw,params,mu);
        itr_in = itr_in + 1;
    end
    
    %% Update d
    eta = eta - (u1 - doR(x,params));
    
    %% End Timer
    itr = itr + 1;

end
