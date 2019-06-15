function w = proj_func(u,v,base_funct,x0,x1)

U = @(x) fevaluation(u,base_funct,x);
V = @(x) fevaluation(v,base_funct,x);

w = ( scprod(U,V,x0,x1)/scprod(U,U,x0,x1)) * u;

end