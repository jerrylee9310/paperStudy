function x1 = run_lsqr
n = 100; 
on = ones(n,1); 
A = spdiags([-2*on 4*on -on],-1:1,n,n);
b = sum(A,2); 
tol = 1e-8; 
maxit = 15;
M1 = spdiags([on/(-2) on],-1:0,n,n); 
M2 = spdiags([4*on -on],0:1,n,n);
x1 = lsqr(@afun,b,tol,maxit,M1,M2);

    function y = afun(x,transp_flag)
       if strcmp(transp_flag,'transp')      % y = A'*x
          y = 4 * x;
          y(1:n-1) = y(1:n-1) - 2 * x(2:n);
          y(2:n) = y(2:n) - x(1:n-1);
       elseif strcmp(transp_flag,'notransp') % y = A*x
          y = 4 * x;
          y(2:n) = y(2:n) - 2 * x(1:n-1);
          y(1:n-1) = y(1:n-1) - x(2:n);
       end
    end
end
