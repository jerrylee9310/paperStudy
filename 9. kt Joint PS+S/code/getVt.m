function [Usinit, tempBasis] = getVt(kdata,L)

    [npe,nfe,nfr] = size(kdata);
    limg = reshape(ifft2(kdata),npe*nfe,nfr);
    [Us,S,Vt] = svds(limg,L);
    Usinit = Us*S;
    tempBasis = Vt';
    
end