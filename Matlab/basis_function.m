function [leg_b,leg_d,dx] = basis_function(N,E2edge,Edge,E2size)

LEG = legendre_coef(N);
LEG_d = poly_deriv(LEG);

legendre_func = @(n,x) poly_eval(n,LEG(n+1,:),x);
legendre_deriv = @(n,x) poly_eval(n-1,LEG_d(n+1,:),x);

affine = @(k,x) 2 / E2size(k) .* (x - (Edge(E2edge(k,1))+Edge(E2edge(k,2)))/2);

leg_b = @(k,n,x) sqrt((2*n+1)/E2size(k)) .*legendre_func(n, affine(k,x));
leg_d = @(k,n,x) 2 / E2size(k) .* sqrt((2*n+1)/E2size(k)) .* legendre_deriv(n, affine(k,x));

dx = @(k1,k2) (E2size(k1)+E2size(k2))/2;

end