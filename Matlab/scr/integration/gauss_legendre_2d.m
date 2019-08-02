function res = gauss_legendre_2d(f,interval,points,weights)

res = 0;

xmin = interval(1);
xmax = interval(2);
ymin = interval(3);
ymax = interval(4);

N = length(points);
points_x = (xmax-xmin)/2 * points + (xmax+xmin)/2;
points_y = (ymax-ymin)/2 * points + (ymax+ymin)/2;
weights_x = weights;
weights_y = weights;

parfor i = 1:N
    for j = 1:N
        res = res + weights_x(i)*weights_y(j)*f(points_x(i),points_y(j));
    end
end

res = res * (xmax-xmin) * (ymax-ymin) / 4;

end