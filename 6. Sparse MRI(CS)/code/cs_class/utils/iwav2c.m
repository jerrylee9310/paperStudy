function res = iwav2c(wdata,widx,filter)
%
%   INPUT :
%       wdata 	: wavelet data
%       widx    : wavelet index
%       filter  : type of filter want to use (ex. 'haar', 'db1'...)
%       
%   OUTPUT : 
%       res : reconstructed data
%

res = waverec2(wdata,widx,filter);
