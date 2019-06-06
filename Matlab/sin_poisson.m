function [f,sol,alpha,mu,dirichlet] = sin_poisson(real)

c = real(1);
d = real(2);

alpha = 1;
mu = 1;

dirichlet = [0 0];

sol = @(x) sin(pi / (d-c) * (x+c));
f = @(x) (1+(pi/(d-c))^2)*sin(pi / (d-c) * (x+c));

end