function [x, u0, u1, u2] = SENSERecon_ALP2(Fhd, eta, u, x, params)
% d: [eta20; eta21; eta22]
% z: [u0;    u1;    u2]  = [Sx;   Rx;   x]
%%
nc = params.nc;
%% Compute SER at current estimate

maxitr = params.maxitr;
maxitr_in = params.maxitr_in;

smap = params.smap;

mu = params.AL.mu;
nu1 = params.AL.nu1;
nu2 = params.AL.nu2;
nuratio = nu2/nu1;

iPpmu = params.AL.iPpmu;
iSpnu2 = params.AL.iSpnu2;
iRpnu2nu1 = params.AL.iRpnu2nu1;

eta20 = eta(:,:,1:nc);
eta22 = eta(:,:,nc+1);
eta21 = eta(:,:,nc+2:end);

u0 = u(:,:,1:nc);
u2 = u(:,:,nc+1);
u1 = u(:,:,nc+2:end);

itr = 1;

%% Inner while loop for minimizing cost corresponding to a given sigma
Sx = smap.*x;
u2 = x;
RW = 0;

while((itr <= maxitr))
    itr_in = 1; % Inner iterations for minimization over z's
    if rem(itr,10) == 0
        disp(['ItrALP2 = ' int2str(itr)]);
    end
    
    while(itr_in <= maxitr_in)
        
        % Compute u0 (Refer revised manuscript)
        rhsz = Fhd + mu*(Sx + eta20);
        u0 = ifft2(iPpmu.*fft2(rhsz));
        
        % Compute u1 (Refer revised manuscript)
        RW = doR(u2, params); % First do Ru2
        Rw = RW + eta21;
        u1 = thresholdRw(Rw,params,nu1*mu);

        % Compute u2 (Refer revised manuscript)
        vmd = u1 - eta21;  % Construct v-dv
        rhsv = doRAdj(vmd, params);                             % doRAdj(v-dv) first
        rhsw = nuratio*(x + eta22) + rhsv;              % Construct the rhs for the system
        u2 = ifft2(iRpnu2nu1.*fft2(rhsw));

        % Compute x (Refer notes)
        rhsx = nu2*(u2 - eta22) + sum(conj(smap) .* (u0 - eta20), 3);
        x = iSpnu2.*rhsx;
        
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
