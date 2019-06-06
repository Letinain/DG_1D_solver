%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is the main script of a 1D DG
% method, with Legendre polynomials and Dirichlet
% conditions.
%
% The equation solved is:
% - mu u'' + alpha u = f
% on the interval [a,b]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear variables;
close all;

% real domain
c = 0; d = 1;
real = [c d];

% [f,sol,alpha,mu,dirichlet_real] = sin_poisson(real);
% [f,sol,alpha,mu,dirichlet_real] = cos_poisson(real);


%%
% simulation domain
a = 0; b = 1;
% simulation = [a b];
simulation = real;
[f,sol,alpha,mu,dirichlet_real] = random_poly_poisson(10,real);
% Dirichlet simulation
ul = 0; ur = rand;
% dirichlet_s = [ul ur];
dirichlet_s = dirichlet_real;



%%
% Arbitrary parameters
beta = 1;
zeta = 0 * [0 1];

%%
% mesh parameters
prec = 20; % number of element
N = 0; % element degree

%%
%%%%%%%%%%%%%%%%%
% mesh generation
%%%%%%%%%%%%%%%%%
[Edge,E2edge,E2size,E2E,normal,K] = mesh_generation(prec,real,simulation,"regular");

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basis function generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%

[leg_b,leg_d,dx] = basis_function(N,E2edge,Edge,E2size);

%%

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

[flux,bound] = flux_matrix(K,N,beta,dirichlet_s,leg_b,dx,leg_d,Edge,E2edge);
% [flux,bound] = flux_matrix_no_bound(K,N,beta,dirichlet_s,leg_b,dx,leg_d,Edge,E2edge);

%%%%%%%%%%%%%
% source term
%%%%%%%%%%%%%

source = source_vector(K,N,f,leg_b,Edge,E2edge);

%%%%%%%%%%%%
% Dirac term
%%%%%%%%%%%%

[Dirac_mat,Dirac_vec] = dirac_term(K,N,real,dirichlet_real,zeta,E2edge,Edge,leg_b);

%%%%%%%%%%
% Assembly
%%%%%%%%%%

A = alpha * mass + mu * (stiff - flux) + Dirac_mat;
F = source + mu* bound + Dirac_vec;

%%
%%%%%%%%%%%%
% resolution
%%%%%%%%%%%%

U = A\F;

%%
%%%%%%%%%%%%%
% total error
%%%%%%%%%%%%%
pt = 10000;
if (~exist('err','var'))
    err = [];
end
err = [err;total_error(N,real,Edge,leg_b,U,pt,sol)];



%%
%%%%%%
% view
%%%%%%
pt_E = 21;
view_simulation(K,N,Edge,E2edge,leg_b,U,pt_E);
hold on;
fplot(sol,real,'r');
view_real(K,N,real,Edge,E2edge,leg_b,U,pt_E);
hold on;
fplot(sol,real,'r')

%%
%%%%%%%%%%%%
% cell error
%%%%%%%%%%%%
pt = 10;
cell_error(K,N,real,Edge,E2edge,leg_b,U,pt,sol);
