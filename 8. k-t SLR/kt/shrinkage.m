function res = shrinkage(A,thres,opt)

switch opt
    case 'sig'
        % singular value shrinkage
        [ny,nx,nt] = size(A);
        
%         [u,s,v] = givefastSVD(reshape(A,ny*nx,nt));
        [u,s,v] = svd(reshape(A,ny*nx,nt),'econ');
        sigv = diag(s) - thres; % threshold
        sigv = (sigv + abs(sigv))/2;
        res = u*diag(sigv)*v';
        res = reshape(res,ny,nx,nt);
        
    case 'mul'
        % multi dimension shrinkage
        mag = sqrt(sos(A));
        mag(mag==0) = 1;
        tmp = max(mag - thres, 0)./mag;
        res = A.*tmp;
        
end