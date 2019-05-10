function smap = getSmap(kdata,lowfreqBand)
% return a coil sensitivity map
% input kdata should have DC point on one of four vertices not center

if nargin == 1
    lowfreqBand = 20;
end

[nx,ny,nc] = size(kdata);

smapMask = zeros(nx,ny);
hanning = hann(lowfreqBand); hanning = hanning.*hanning.';

smapMask(floor((nx-lowfreqBand)/2)+1:floor((nx+lowfreqBand)/2),floor((nx-lowfreqBand)/2)+1:floor((nx+lowfreqBand)/2)) = ...
    hanning;
smapMask = fftshift(smapMask);
lowImg = ifft2(smapMask.*kdata);
lowsosImg = sqrt(sum(abs(lowImg).^2,3)); 
smap = lowImg./lowsosImg;