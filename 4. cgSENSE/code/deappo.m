function gx=deappo(b,wid,im_size)
% function deap_func=deappo(kern_const,kern_wid,im_size,eps)
%% deappodization function for convolution gridding

% only for kaiser-bessel kernel

x = -im_size/2:im_size/2-1;
cx = sinc(sqrt((pi*wid*x/im_size).^2-b^2)/pi)/im_size;
gx = cx'*cx;

% om = [0:kern_wid/2]/(kern_wid/2);
% kaiser_kern = besseli(0,kern_const*sqrt(1-om.*om));
% kaiser_kern = kaiser_kern./kaiser_kern(1);
% k_kern = [kaiser_kern(end:-1:2) kaiser_kern];
% 
% deap_func = zeros(im_size);
% deap_func(round((im_size-kern_wid)/2 +1):round((im_size+kern_wid)/2),round((im_size-kern_wid)/2 +1):round((im_size+kern_wid)/2)) = k_kern'*k_kern;
% 
% deap_func = abs(fft2c(deap_func)) + eps;