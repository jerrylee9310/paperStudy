clear all; clc; close all;
addpath(strcat(pwd,'/utils')); 

img = phantom(100);

sens = gausswin(160) * gausswin(160)';

% coil sensitivity
cSens = zeros(200,200,4);
cSens(1:160,21:180,1) = sens;
cSens(21:180,41:end,2) = sens;
cSens(41:end,21:180,3) = sens;
cSens(21:180,1:160,4) = sens;
cSens = cSens(51:150,51:150,:);

% image data
imgc = img.*cSens;

kdata = i






tt = imgc;
for i = 1:4
    subplot(2,2,i); imshow(abs(tt(:,:,i)),[]);
end