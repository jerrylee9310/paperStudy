function res = mtimes(A,x)

if A.adjoint == 0
    switch A.type
        case 'W'
            ...
            
        case 'D'
            % Dx
            Dx = imfilter(x,[-1 1],'circular');

            % Dy
            Dy = imfilter(x,[-1 1]','circular'); 

            % Dt
            Dt = imfilter(x,reshape([-1 1],1,1,2),'circular');
            
            res = cat(4,Dx,Dy,Dt);
            
        case 'WD'
            ...
                
    end
    
    ...
else
    switch A.type
        case 'W'
            ...
                
        case 'D'
            Dx = x(:,:,:,1);
            Dy = x(:,:,:,2);
            Dt = x(:,:,:,3);
            
            % iDx
            iDx = cat(2,Dx(:,end,:),Dx);
            iDx = imfilter(iDx,[1 -1], 'circular');
            iDx = iDx(:,1:end-1,:);
            

            % iDy
            iDy = cat(1,Dy(end,:,:),Dy);
            iDy = imfilter(iDy,[1 -1]', 'circular');
            iDy = iDy(1:end-1,:,:);
           
            % iDt     
            iDt = cat(3,Dt(:,:,end),Dt);
            iDt = imfilter(iDt,reshape([1 -1],1,1,2), 'circular');
            iDt = iDt(:,:,1:end-1);
            
            res = iDx + iDy + iDt;
            
        case 'WD'
            ...
                
    end
end