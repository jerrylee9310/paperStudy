function res = FTS(x,mask,opt)

    if nargin == 2
        opt = 'forward';
    end

    [np,nf,nt] = size(mask);
    
    switch opt
        case 'forward'
            x = reshape(x,size(mask));
            uk = (1/sqrt(np*nf))*fft2(x);
            res = uk(mask>0);
        case 'adjoint'
            mk = zeros(np,nf,nt);
            mk(mask>0) = x;
            res = sqrt(np*nf)*ifft2(mk);
    end
    
    res = res(:);