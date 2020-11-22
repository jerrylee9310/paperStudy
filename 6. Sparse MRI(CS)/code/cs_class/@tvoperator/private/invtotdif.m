function res = invtotdif(y)

res = adjDx(y(:,:,1)) + adjDy(y(:,:,2));

return;

function res = adjDy(x)
res = x(:,[1,1:end-1]) - x;
res(:,1) = -x(:,1);
res(:,end) = x(:,end-1);

function res = adjDx(x)
res = x([1,1:end-1],:) - x;
res(1,:) = -x(1,:);
res(end,:) = x(end-1,:);


