function y = dirichlet_flux_interior(beta,leg_b,leg_d,dx,k1,i,j,x,n)

y = (-beta/dx(k1,k1) * leg_b(k1,j,x) + n * leg_d(k1,j,x) )* leg_b(k1,i,x);

end