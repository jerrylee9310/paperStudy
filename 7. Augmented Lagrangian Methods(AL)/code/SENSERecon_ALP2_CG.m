function [x] = SENSERecon_ALP2_TEST(Fhd, eta, u, x, params)
% d: [eta20; eta21; eta22]
% z: [u0;    u1;    u2]  = [Sx;   Rx;   x]
%% Compute SER at current estimate

maxitr = params.maxitr;
maxitr_in = params.maxitr_in;

smap = params.smap;
mask = params.mask;

mu = params.AL.mu;
nu1 = params.AL.nu1;
nu2 = params.AL.nu2;
nuratio = nu2/nu1;

iPpmu = params.AL.iPpmu;
iSpnu2 = params.AL.iSpnu2;
iRpnu2nu1 = params.AL.iRpnu2nu1;

eta20 = eta(:,:,1:12);
eta21 = eta(:,:,13:21);
eta22 = eta(:,:,22);

u0 = u(:,:,1:12);
u1 = u(:,:,13:21);
u2 = u(:,:,22);

itr = 1;

%% Inner while loop for minimizing cost corresponding to a given sigma
Sx = smap.*x;
u2 = x;
RW = 0;

while((itr <= maxitr))
    itr_in = 1; % Inner iterations for minimization over z's
%     if rem(itr,10) == 0
        disp(['ItrALP2 = ' int2str(itr)]);
%     end
    
    while(itr_in <= maxitr_in)
        
        % Compute u0 (Refer revised manuscript)
        rhsz = Fhd + mu*(Sx + eta20);
        
%         u0 = ifft2(iPpmu.*fft2(rhsz));
        
        r = rhsz - (ifft2(mask.*fft2(u0)) + mu*u0);
        p = r;
        err0 = r(:)'*r(:);

        %% main loop
        for i = 1:100
            aNu = r(:)'*r(:);
            Ap = ifft2(mask.*fft2(p)) + mu*p;
            aDe = abs(p(:)'*Ap(:)); % p'*A*p
            alpha = aNu / aDe;

            u0 = u0 + alpha*p;
            r_nxt = r - alpha*Ap; % r - alpha*A*p

            beta = (r_nxt(:)'*r_nxt(:)) / (r(:)'*r(:));

            p = r_nxt + beta*p;

            err = r_nxt(:)'*r_nxt(:);
            if err/err0 < 1e-100
                break
            end

            r = r_nxt;
        end
        
        % Compute u1 (Refer revised manuscript)
        RW = doR(u2, params); % First do Ru2
        Rw = RW + eta21;
        u1 = thresholdRw(Rw,params,nu1*mu);

        % Compute u2 (Refer revised manuscript)
        vmd = u1 - eta21;  % Construct v-dv
        rhsv = doRAdj(vmd, params);                             % doRAdj(v-dv) first
        rhsw = nuratio*(x + eta22) + rhsv;              % Construct the rhs for the system
        
%         u2 = ifft2(iRpnu2nu1.*fft2(rhsw));
        
        r = rhsw - (doRAdj(doR(u2,params),params) + nuratio*u2);
        p = r;
        err0 = r(:)'*r(:);

        %% main loop
        for i = 1:100
            aNu = r(:)'*r(:);
            Ap = doRAdj(doR(p,params),params) + nuratio*p;
            aDe = abs(p(:)'*Ap(:)); % p'*A*p
            alpha = aNu / aDe;

            u2 = u2 + alpha*p;
            r_nxt = r - alpha*Ap; % r - alpha*A*p

            beta = (r_nxt(:)'*r_nxt(:)) / (r(:)'*r(:));

            p = r_nxt + beta*p;

            err = r_nxt(:)'*r_nxt(:);
            if err/err0 < 1e-30
                break
            end

            r = r_nxt;
        end
        

        % Compute x (Refer notes)
        rhsx = nu2*(u2 - eta22) + sum(conj(smap) .* (u0 - eta20), 3);
        
%         x = iSpnu2.*rhsx;
        
        r = rhsx - ( sum(conj(smap).*smap.*x,3) + nu2*x );
        p = r;
        err0 = r(:)'*r(:);

        %% main loop
        for i = 1:100
            aNu = r(:)'*r(:);
            Ap = sum(conj(smap).*smap.*p,3) + nu2*p;
            aDe = abs(p(:)'*Ap(:)); % p'*A*p
            alpha = aNu / aDe;

            x = x + alpha*p;
            r_nxt = r - alpha*Ap; % r - alpha*A*p

            beta = (r_nxt(:)'*r_nxt(:)) / (r(:)'*r(:));

            p = r_nxt + beta*p;

            err = r_nxt(:)'*r_nxt(:);
            if err/err0 < 1e-150
                break
            end

            r = r_nxt;
        end
        
        
        % Book Keeping
        Sx = smap .* x;
        RW = doR(u2, params);

        itr_in = itr_in + 1;
    end
    
    %% Update d
    eta20 = eta20 - (u0 - Sx);
    eta21 = eta21 - (u1 - RW);
    eta22 = eta22 - (u2 - x);
    
    %% End Timer
    itr = itr + 1;

end
