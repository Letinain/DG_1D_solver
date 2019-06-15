function n = fnorm(u,basis_funct,x0,x1)

f = @(x) fevaluation(u,basis_funct,x);

n = sqrt(scprod(f,f,x0,x1));

end