function mask = genmask(dsize,maskparams)
%% generate sampling mask

if length(dsize) == 2
    % 2D mask
    
elseif length(dsize) == 3
    % 2D+time mask
%     [ny,nx,ny] = dsize;
    np = dsize(1); nf = dsize(2); nt = dsize(3);
    switch maskparams.type
        case 'c'
            lwb = maskparams.lowfreq;
            rfac = maskparams.rfac;
            
            temmask = zeros(np, nt);
            Nsam = round(((nf*rfac - lwb)*np)/(np-lwb));
            nav_ind = np/2-lwb/2+1:np/2+lwb/2;
            for ii = 1:np
                if (find(ii == nav_ind))
                    temmask(ii, :) = ones(1, nt);
                else
                    temp = zeros(1, nt);
                    ind = randperm(nt);
                    temp(ind(1:Nsam)) = 1;
                    temmask(ii, :) = temp;
                end
            end
            temmask = logical(temmask);
            mask = permute(reshape(repmat(temmask, 1, nf), np, nt, nf), [1, 3, 2]);
                
        case 'r'
            nspk = maskparams.nspoke;
            ny = np; nx = nf;
            for t = 1:nt
                y = -ny/2:1:ny/2-0.5;
                x = -nx/2:1:nx/2-0.5;
                
                % slight different radial angle between PE (Dt needs
                % incoherence sampling along time direction)
                dtht = (pi/nspk)*rand;
                
                mm = [];
                for ang = (0:pi/nspk:pi-1e-3)
                    map = complex(y*cos(ang+dtht),x*sin(ang+dtht));
                    mm = [mm; map];
                end
                
                % round the collected data to the nearest cartesian
                % location
                kcart = round(mm+(0.5+0.5i));
                % plot(kcart,'*');title('k locations after nearest neighbor interpolation: Center (0,0)');
                
                % shift the cartesian locations 
                kloc1 = round(kcart)+((ny/2+1)+(nx/2+1)*1i);
                kloc1real = real(kloc1); kloc1real = kloc1real - ny*(kloc1real>ny);
                kloc1imag = imag(kloc1); kloc1imag = kloc1imag - nx*(kloc1imag>nx);
                kloc1real = kloc1real + ny*(kloc1real<1);
                kloc1imag = kloc1imag + nx*(kloc1imag<1);
                kloc1 = kloc1real + 1i*kloc1imag;
                
                % sampling pattern
                for ii=1:size(kloc1,1)
                    for j=1:size(kloc1,2)

                    mask(real(kloc1(ii,j)),imag(kloc1(ii,j)),t) = 1;

                    end
                end
            end
        
    end
end
