function res = FTT(x,size,opt)

    if nargin == 2
        opt = 'forward';
    end

    np = size(1); nf = size(2); nt = size(3);
    
    switch opt
        case 'forward'
            td = reshape(x,[np nf nt]);
            res = (1/sqrt(nt))*fft(td,[],3);
        case 'adjoint'
            fd = reshape(x,[np nf nt]);
            res = sqrt(nt)*ifft(fd,[],3);
    end
    res = res(:);