function [f,sol,alpha,mu,neumann] = cos_neumann(real)

c = real(1);
d = real(2);

alpha = 1;
mu = 1;

k = 5;

neumann = [0 0];

sol = @(x) cos(k * pi * (x-c) / (d-c));
f = @(x) (alpha + mu * k * pi * k * pi) * cos(k * pi * (x-c) / (d-c));

end