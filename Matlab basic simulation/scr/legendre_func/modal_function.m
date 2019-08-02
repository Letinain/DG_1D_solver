function y = modal_function(f,n,Edge,sol,x)
y = zeros(size(x));

[i,j] = size(x);

for l=1:i
    for m=1:j
        k = 1;
        while (x(l,m)>Edge(k+1))
            k = k+1;
        end
        
        for i=0:n
            y(l,m) = y(l,m) + sol((k-1)*(n+1)+i+1) * f(k,i,x(l,m));
        end
    end
end
end