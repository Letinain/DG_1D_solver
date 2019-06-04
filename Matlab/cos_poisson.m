function [f,sol,alpha,mu,dirichlet] = cos_poisson(real)

c = real(1);
d = real(2);

alpha = 1/2;
mu = 2*((d-c)/pi)^2;

dirichlet = [1 0];

sol = @(x) cos(pi/2 * (x-c) / (d-c));
f = @(x) cos(pi/2 * (x-c) / (d-c));

end