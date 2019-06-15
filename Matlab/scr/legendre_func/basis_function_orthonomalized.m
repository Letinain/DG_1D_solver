function [leg_b,leg_d,dx] = basis_function_orthonomalized(N,E2edge,Edge,E2size,real)

[leg_b,leg_d,~] = basis_function(N,E2edge,Edge,E2size);

k = 1;
coeffs1 = speye(N);
basis_funct = @(n,x) leg_b(k,n-1,x);
x0 = real(1); x1 = Edge(2);

for ii= 1:N
    coeffs1(:,ii) = coeffs1(:,ii)/fnorm(coeffs1(:,ii),basis_funct,x0,x1);
    for jj = ii+1:N
        coeffs1(:,jj) = coeffs1(:,jj) - proj_func(coeffs1(:,ii),coeffs1(:,jj),basis_funct,x0,x1);
    end
end

k = 2;
coeffs2 = speye(N);
basis_funct = @(n,x) leg_b(k,n-1,x);
x0 = real(1); x1 = Edge(2);

for ii= 1:N
    coeffs2(:,ii) = coeffs2(:,ii)/fnorm(coeffs2(:,ii),basis_funct,x0,x1);
    for jj = ii+1:N
        coeffs2(:,jj) = coeffs2(:,jj) - proj_func(coeffs2(:,ii),coeffs2(:,jj),basis_funct,x0,x1);
    end
end





%%

len = Edge(2:end) - Edge(1:end-1);
len(1) = Edge(2)-real(1);
len(end) = real(2)-Edge(end-1);

dx = @(k1,k2) (len(k1)+len(k2))/2;



end
