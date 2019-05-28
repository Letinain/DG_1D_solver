function y = legendre_derivative(a,b,n,x)
    h = abs(b-a);
    moy = (a+b) / 2;
    x_ref = 2 / h * (x - moy);
    y = 0;
    if (mod(n,2)==0)
        for i=0:n/2-1
            y = y + (4*i+3) * legendreP(2*i+1,x_ref);
        end
    else
        for i = 0:(n-1)/2
            y = y + (4*i+1) * legendreP(2*i,x_ref);
        end
    end

    y = y * sqrt((2*n+1)/2) * (2/h)^(3/2);
    
end