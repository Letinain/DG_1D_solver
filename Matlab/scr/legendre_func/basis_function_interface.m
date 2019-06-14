function [leg_b,leg_d,dx] = basis_function_interface(N,E2edge,Edge,E2size,real)

[leg_b,leg_d,~] = basis_function(N,E2edge,Edge,E2size);

len = Edge(2:end) - Edge(1:end-1);
len(1) = Edge(2)-real(1);
len(end) = real(2)-Edge(end-1);

dx = @(k1,k2) (len(k1)+len(k2))/2;

end