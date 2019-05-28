function y = legendre_basis(a,b,n,x)
    h = abs(b-a);
    moy = (a+b) / 2;
    y = sqrt((2*n+1)/h) * legendreP(n, 2 / h * (x - moy) );
end