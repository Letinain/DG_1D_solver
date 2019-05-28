function y = poly_eval(n,A,x)

y = A(1);
X = x;
for i=1:n
    y = y + A(1,i+1)*X;
    X = X.*x;
end


end