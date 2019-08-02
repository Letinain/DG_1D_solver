function [leg_b,leg_d,dx] = basis_function(N,E2edge,Edge,E2size)

LEG = zeros(N+1);
for i = 0:N
    LEG(i+1,end-i:end) = LegendrePoly(i);
end

LEG_d = zeros(N+1,N);
for i = 0:N
    tmp = polyder(LEG(i+1,:));
    LEG_d(i+1,end-length(tmp)+1:end) = tmp;
end

legendre_func = @(n,x) polyval(LEG(n+1,:),x);
legendre_deriv = @(n,x) polyval(LEG_d(n+1,:),x);

affine = @(k,x) 2 / E2size(k) .* (x - (Edge(E2edge(k,1))+Edge(E2edge(k,2)))/2);

leg_b = @(k,n,x) sqrt((2*n+1)/E2size(k)) .*legendre_func(n, affine(k,x));
leg_d = @(k,n,x) 2 / E2size(k) .* sqrt((2*n+1)/E2size(k)) .* legendre_deriv(n, affine(k,x));

dx = @(k1,k2) (E2size(k1)+E2size(k2))/2;

end