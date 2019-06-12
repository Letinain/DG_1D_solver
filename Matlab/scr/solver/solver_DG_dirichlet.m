function U = solver_DG_dirichlet(K,N,alpha,mu,beta,dirichlet_s,leg_b,dx,leg_d,Edge,E2edge,E2size,E2E,Edge2bound,normal,f)

%%%%%%%%%%%%%
% mass matrix
%%%%%%%%%%%%%

mass = mass_matrix(K,N);

%%%%%%%%%%%%%%%%%%
% stiffness matrix
%%%%%%%%%%%%%%%%%%

stiff = stiff_matrix(K,N,E2size);

%%%%%%%%%%%%%
% flux matrix
%%%%%%%%%%%%%

[flux,bound] = flux_matrix_dirichlet(K,N,beta,dirichlet_s,leg_b,dx,leg_d,Edge,E2edge,E2E,Edge2bound,normal);

%%%%%%%%%%%%%
% source term
%%%%%%%%%%%%%

source = source_vector(K,N,f,leg_b,Edge,E2edge);

%%%%%%%%%%
% Assembly
%%%%%%%%%%

A = alpha * mass + mu * (stiff - flux);
F = source + mu* bound;


%%%%%%%%%%%%
% resolution
%%%%%%%%%%%%

U = A\F;


end