function res = getPsi(A,params)
    switch params.prior.type
        case 'D'
%             tmp1 = norm(xtogam(A(:,:,:,1)),'fro').^2;
%             tmp2 = norm(reshape(A(:,:,:,2),128*128,50),'fro').^2;
%             tmp3 = norm(reshape(A(:,:,:,3),128*128,50),'fro').^2;
%             res = sqrt(tmp1+tmp2+tmp3);
            tmp = sqrt(abs(A(:,:,:,1)).^2 + abs(A(:,:,:,2)).^2 + abs(A(:,:,:,3)).^2);
            res = sum(tmp(:));
    end
    