function res = mtimes(A,x)

if A.adjoint == 1
    % inverse Fourier transform on k-space
    res = ifft2(A.mask.*x);
    
else
    % Fourier transform on image
    res = A.mask.*fft2(x);
    
end
