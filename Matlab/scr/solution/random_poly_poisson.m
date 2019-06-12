function [f,sol,alpha,mu,dirichlet] = random_poly_poisson(r,real)

P = 10*rand(1,r+1)-5;
P_d = poly_deriv(P);

c = real(1);
d = real(2);

alpha = 0.5 + rand;
mu = 0.5 + rand;

dirichlet = [0 0];

P_dd = poly_deriv(P_d);
sol = @(x) x.*(x-1).*poly_eval(r,P,x);

f = @(x) (alpha*sol(x) - mu * (x.*(x-1).*poly_eval(r-2,P_dd,x) + (4*x-2).*poly_eval(r-1,P_d,x) + 2*poly_eval(r,P,x))).*(x<=d) .* (x>=c);
% f = @(x) f_mid(x) - f_mid(c).*(x<c) - f_mid(d).*(x>d);
end