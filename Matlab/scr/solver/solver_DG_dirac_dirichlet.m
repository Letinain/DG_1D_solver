function U = solver_DG_dirac_dirichlet(K,N,alpha,mu,beta,dirichlet_s,leg_b,dx,leg_d,Edge,E2edge,E2size,f,real,dirichlet,zeta)

% mass matrix
mass = mass_matrix(K,N);

% stiffness matrix
stiff = stiff_matrix(K,N,E2size);

% flux matrix
[flux,bound] = flux_matrix_dirichlet(K,N,beta,dirichlet_s,leg_b,dx,leg_d,Edge,E2edge);

% source term
source = source_vector(K,N,f,leg_b,Edge,E2edge);

% dirac term
[dirac_mat,dirac_vec] = dirac_term(K,N,real,dirichlet,zeta,E2edge,Edge,leg_b);

% Assembly
A = alpha * mass + mu * (stiff - flux) + dirac_mat;
F = source + mu* bound + dirac_vec;

% resolution
U = A\F;


end