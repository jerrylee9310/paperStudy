function [wdata,widx] = wav2c(x,level,filter)
%
%   INPUT :
%       x   	: input data
%       level   : wavelet depth level
%       filter  : type of filter want to use (ex. 'haar', 'db1'...)
%       
%   OUTPUT : 
%       c : wavelet data
%       s : level-index data
%

[wdata,widx] = wavedec2(x,level,filter);
wdata = wdata.';
