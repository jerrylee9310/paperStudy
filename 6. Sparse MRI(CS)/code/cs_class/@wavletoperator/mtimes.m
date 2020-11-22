function res = mtimes(W,x)

if W.adjoint == 1
    % Phi'*wav = img
    % x = wavelet domain(vector)
    res = iwav2c(x,W.idx,W.filter);
else
    % Phi*img = wav
    % img = matrix
    [res,idx] = wav2c(x,W.level,W.filter);
end
