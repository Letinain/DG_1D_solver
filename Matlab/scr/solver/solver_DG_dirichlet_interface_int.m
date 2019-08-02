function [U,A,F,mass,stiff,flux,source,bound] = solver_DG_dirichlet_interface_int(K,N,real,alpha,mu,beta,f,dirichlet,leg_b,leg_d,dx,Edge,E2edge,E2E,E2size,E2bound,normal)

% gauss-legendre integration points/weight
[x,w]=lgwt(N+1,-1,1);
gauss = [x';w'];

% mass matrix
mass = mass_matrix_interface_int(K,N,real,leg_b,Edge,E2edge,E2bound,gauss);

% stiffness matrix

stiff = stiff_matrix_interface_int(K,N,real,leg_d,Edge,E2edge,E2size,E2bound,gauss);

% flux matrix
[flux,bound] = flux_matrix_interface(K,N,real,beta,dirichlet,leg_b,dx,leg_d,Edge,E2edge,E2E,E2bound,normal);

% source term
source = source_vector_interface_int(K,N,f,leg_b,real,Edge,E2edge,E2bound,gauss);

% Assembly
A = alpha*mass + mu*(stiff - flux);
F = source + mu*bound;

% resolution
U = A\F;

end