function err = total_error_int(N,domain,Edge,leg_b,U,sol)

c = domain(1);
d = domain(2);

[x,w]=lgwt(N+1,c,d);

approx = @(x) modal_function(leg_b,N,Edge,U,x);
err = zeros(3,2);
err(1,1) = gauss_legendre_1d(@(x) abs(approx(x)-sol(x)),x,w);
err(1,2) = 100 * err(1,1)/gauss_legendre_1d(@(x) abs(sol(x)),x,w);
err(2,1) = sqrt(gauss_legendre_1d(@(x) (approx(x)-sol(x)).^2,x,w));
err(2,2) = 100 * err(2,1)/sqrt(gauss_legendre_1d(@(x) sol(x).^2,x,w));
err(3,1) = max(abs(approx(x)-sol(x)));
err(3,2) = 100 * err(3,1)/max(abs(sol(x)));
% fprintf('error L1: %e (%2.6f %%) \nerror L2: %e (%2.6f %%)\nerror max: %e (%2.6f %%)\n',err(1,1),err(1,2),err(2,1),err(2,2),err(3,1),err(3,2)); 

err = err(:,1)';

end