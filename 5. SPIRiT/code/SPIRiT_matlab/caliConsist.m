function res = caliConsist(kernel,data)

ncoil = 4;
ns = 256;
res = zeros(size(data));

for nc=1:ncoil
    convker = kernel(:,:,:,nc);
    for nx=1:ns
        for ny = 1:ns
            res(nx,ny) = sum(sum(convker.*data(nx:nx+4,ny:ny+4,:)));
end