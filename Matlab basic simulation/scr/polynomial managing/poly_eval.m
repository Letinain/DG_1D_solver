function y = poly_eval(A,x)

y = A(end);
N = length(A);
for i = N-1:-1:1
    y = A(i) + x .* y;
end


% y = A(1) * ones(size(x));
% X = x;
% for i=1:n
%     y = y + A(1,i+1)*X;
%     X = X.*x;
% end


end