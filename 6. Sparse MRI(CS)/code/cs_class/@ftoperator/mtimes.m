function res = mtimes(FT,x)

if FT.adjoint == 1
    % FT'(k) = img
    % k : k-space data
    % FT' : partial inverse fourier transform (undersampling ifft)
    % img : image
    
    tmp = reshape(x,FT.imSize);
    res = tmp.*FT.mask;
    res = ifft2c(res);

else
    % FT(img) = k
    % img : image
    % FT : partial inverse fourier transform (undersampling fft)
    % k : undersampled k-space
    
    tmp = reshape(x,FT.imSize);
    res = fft2c(tmp);
    res = res.*FT.mask;
end

