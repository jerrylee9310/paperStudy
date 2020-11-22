function res = wavletoperator(img, filter, level)

res.adjoint = 0;
res.level = level;
res.filter = filter;
[nothing, res.idx] = wav2c(img, level, filter);

res = class(res,'wavletoperator');