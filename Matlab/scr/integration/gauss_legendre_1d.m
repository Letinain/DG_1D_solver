function res = gauss_legendre_1d(f,points,weights)

res = 0;
for i=1:length(points)
    res = res + weights(i)*f(points(i));
end

end