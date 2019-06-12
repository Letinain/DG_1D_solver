function [f,sol,alpha,mu,dirichlet] = sin_poisson(k,real)

c = real(1);
d = real(2);

alpha = 1;
mu = 1;

dirichlet = [0 0];

sol = @(x) sin(k * pi / (d-c) * (x-c));
f = @(x) (1+ k^2 * (pi/(d-c))^2)*sin(k * pi / (d-c) * (x-c));

end