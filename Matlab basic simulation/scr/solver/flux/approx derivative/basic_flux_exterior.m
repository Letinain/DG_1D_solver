function y = basic_flux_exterior(beta,leg_b,leg_d,dx,k1,k2,i,j,x,n)

y = (  beta / dx(k1,k2) * leg_b(k2,j,x) + n * leg_d(k2,j,x) / 2) * leg_b(k1,i,x);

end