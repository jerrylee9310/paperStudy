function res = wavshow(wdata,widx,filter)

lvl = size(widx,1) - 2;

for i = 1:lvl
    [H{i}, V{i}, D{i}] = detcoef2('all',wdata,widx,i);
end

A = appcoef2(wdata,widx,filter,lvl);
res = A;

for i = lvl:-1:1
    res = [res H{i}; V{i} D{i}];
end

% figure(); imshow(abs(res),[]);