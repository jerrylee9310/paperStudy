function res = sparsity_prior_operator(params)

res.adjoint = 0;
res.type = params.prior.type;

res = class(res,'sparsity_prior_operator');