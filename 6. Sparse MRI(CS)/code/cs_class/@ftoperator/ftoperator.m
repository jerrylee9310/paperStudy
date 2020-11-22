function res = ftoperator(mask)

res.adjoint = 0; % transpose flag
res.mask = mask;
res.imSize = size(mask);

res = class(res,'ftoperator');