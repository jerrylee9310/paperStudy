function dc_idx = findDC(kdata)
%%%%%%
%-----INPUT
%       Cartesian data
%
%%%%%%
if ndims(kdata) == 2
    [nx,ny] = size(kdata);
    nc = 1;
else
    [nx,ny,nc] = size(kdata);
end

dc_idx = zeros(1,2,nc);
for c = 1:nc
    data = kdata(:,:,c);
    [~,max_idx] = max(data(:));
    [x_idx,y_idx] = ind2sub(size(data),max_idx);
    dc_idx(:,:,c) = [x_idx,y_idx];
end
