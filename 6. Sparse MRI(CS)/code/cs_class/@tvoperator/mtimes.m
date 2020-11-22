function res = mtimes(TV,x)


if TV.adjoint
	res = invtotdif(x);

else
	res = totdif(x);
    
end