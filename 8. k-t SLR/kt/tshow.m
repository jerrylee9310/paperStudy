function tshow(data)

[nx,ny,nt] = size(data);


figure(),
for t = 1:nt
    imshow(abs(data(:,:,t)),[]); pause(0.01);
end
