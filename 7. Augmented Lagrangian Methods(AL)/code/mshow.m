function mshow(data,row,col,slice)

if nargin == 1
	row = 1;
    col = 1;
    slice = 1;
elseif nargin == 3
        slice = row*col;
end

figure(),
for i = 1:slice
    subplot(row,col,i); imshow(abs(data(:,:,i)),[]);
end
