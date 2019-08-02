function [f,sol,alpha,mu,dirichlet] = sin_poisson(k,real)

c = real(1);
d = real(2);

alpha = 0;
mu = 1;

dirichlet = [1 1];

sol = @(x) sin(k * pi / (d-c) * (x-c))+1;
f = @(x) (k^2 * (pi/(d-c))^2)*sin(k * pi / (d-c) * (x-c));

end