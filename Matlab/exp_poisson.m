function [f,sol,alpha,mu,dirichlet] = exp_poisson(real)

c = real(1);
d = real(2);

alpha = 1;
mu = 1;

dirichlet = [0 0];

sol = @(x) (x-1).*(exp(-x)-1);
f = @(x) (sol(x) - (x-3).*exp(-x)) .* (x<=d) .* (x>=c) - (x>d) .* (sol(d) - (d-3).*exp(-d))...
    -(x<c) .* (sol(c) - (c-3).*exp(-c));

end