function res = totdif(image)

[sx,sy] = size(image);

Dx = image([2:end,end],:) - image;
Dy = image(:,[2:end,end]) - image;

res = cat(3,Dx,Dy);