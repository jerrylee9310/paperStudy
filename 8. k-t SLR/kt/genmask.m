function mask = genmask(dsize,maskparams)

if length(dsize) == 2
    % 2D mask
    
elseif length(dsize) == 3
    % 2D+time mask
%     [ny,nx,ny] = dsize;
    ny = dsize(1); nx = dsize(2); nt = dsize(3);
    switch maskparams.type
        case 'c'
            ...
                
        case 'r'
            nspk = maskparams.nspoke;
            
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
