function U = solver_DG_neumann(K,N,alpha,mu,beta,neumann,leg_b,dx,leg_d,Edge,E2edge,E2size,f)

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

[flux,bound] = flux_matrix_Neuman(K,N,beta,neumann,leg_b,dx,leg_d,Edge,E2edge);

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