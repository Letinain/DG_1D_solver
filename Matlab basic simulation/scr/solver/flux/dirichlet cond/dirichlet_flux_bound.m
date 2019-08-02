function y = dirichlet_flux_bound(beta,leg_b,dx,k,i,x,diri)

y = beta / dx(k,k) * diri * leg_b(k,i,x);

end