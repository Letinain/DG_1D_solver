function [f,sol,alpha,mu,dirichlet] = sin_poisson(real)

c = real(1);
d = real(2);

alpha = 1/2;
mu = ((d-c)/pi)^2 / 2;

dirichlet = [0 0];

sol = @(x) sin(pi / (d-c) * (x+c));
f = @(x) sin(pi / (d-c) * (x+c));

end