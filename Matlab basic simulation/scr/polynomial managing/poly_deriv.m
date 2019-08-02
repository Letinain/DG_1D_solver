function A = poly_deriv(A)
[n,m] = size(A);
for i = 0:m-2;
    A(:,i+1) = (i+1)*A(:,i+2);
end
A(:,end) = 0;
end