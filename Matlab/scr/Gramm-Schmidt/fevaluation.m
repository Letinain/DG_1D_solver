function y = fevaluation(u,base_funct,x)
    y = 0;
    n = length(u);
    for i = 1:n
        y = y + u(i) * base_funct(i,x);
    end
end