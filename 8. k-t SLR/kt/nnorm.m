function res = nnorm(A)
    
    [ny,nx,nt] = size(A);
    if length(size(A)) == 3
        A = reshape(A,ny*nx,nt);
    end
    
    s = svd(A);
    sv = diag(s);
    res = sum(abs(sv(:)));