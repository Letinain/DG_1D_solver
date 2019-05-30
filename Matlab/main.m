%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is the main script of a 1D DG
% method, with Legendre polynomials and Dirichlet
% conditions.
%
% The equation solved is:
% - alpha u'' + mu u = f
% on the interval [a,b]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables;
close all;

% simulation domain
a = -100; b = 101;
simulation = [a b];

% real domain
c = 0; d = 1;
real = [c d];

% Problem parameters
alpha = 1; mu = 1;
%mu = 1/2 * ((d-c)/pi)^2;
%mu = 2/pi^2;

% Arbitrary parameters
beta = 1;
zeta = [-1000 -1000];

% Dirichlet
ul = 0; ur = 0;
dirichlet = [ul ur];

% solution
%sol = @(x) (x-1).*(exp(-x)-1);
%sol = @(x) sin(pi*x);
%sol = @(x) cos(pi/2*x);

r = 10;
P = 10*rand(1,r+1)-5;
P_d = poly_deriv(P);
P_dd = poly_deriv(P_d);
sol = @(x) x.*(x-1).*poly_eval(r,P,x);



% source
%f = @(x) (sol(x) - (x-3).*exp(-x)) .* (x<=d) .* (x>=c) - (x>d) .* (sol(d) - (d-3).*exp(-d))...
%    -(x<c) .* (sol(c) - (c-3).*exp(-c));
%f = @(x) sin(pi / (d-c) * (x+c));
%f = @(x) cos(pi/2 * (x+a) / (b-a));
f = @(x) (mu*sol(x) - alpha * (x.*(x-1).*poly_eval(r-2,P_dd,x) + (4*x-2).*poly_eval(r-1,P_d,x) + 2*poly_eval(r,P,x))).*(x<=d) .* (x>=c);

% mesh parameters
K = 800; % number of element
N = 3; % element degree

%%
%%%%%%%%%%%%%%%%%
% mesh generation
%%%%%%%%%%%%%%%%%
[Edge,E2edge,E2size,E2E,normal] = mesh_generation(K,simulation,"regular");

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

[flux,bound] = flux_matrix(K,N,beta,dirichlet,leg_b,dx,leg_d,Edge,E2edge);

%%%%%%%%%%%%%
% source term
%%%%%%%%%%%%%

source = source_vector(K,N,f,leg_b,Edge,E2edge);

%%%%%%%%%%%%
% Dirac term
%%%%%%%%%%%%

Dirac = dirac_vector(K,N,real,zeta,E2edge,Edge,leg_b);

%%%%%%%%%%
% Assembly
%%%%%%%%%%

A = alpha * mass + mu * (stiff - flux) + Dirac;
F = source + mu* bound;

%%%%%%%%%%%%
% resolution
%%%%%%%%%%%%

U = A\F;

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
%%%%%%%%%%%%%
% total error
%%%%%%%%%%%%%
pt = 10000;
total_error(N,real,Edge,leg_b,U,pt,sol);

%%
%%%%%%%%%%%%
% cell error
%%%%%%%%%%%%
pt = 10;
cell_error(K,N,real,Edge,E2edge,leg_b,U,pt,sol);
