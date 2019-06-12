addpath(genpath(pwd));

% real domain
c = -1000; d = 1000;
real = [c d];
l = abs(d-c);
center = (c+d)/2;

% simulation domain
a = -1; b = 2;
simulation = [a b];
dirichlet_s = [0 0];

% n = 1;
% [f,sol,alpha,mu,dirichlet] = random_poly_poisson(n,real);
[f,sol,alpha,mu,dirichlet] = sin_poisson(5,real);
% n = 1;
% sol = @(x) x.^n;
% f = @(x) sol(x)-n*(n-1)*x.^(n-2);
% alpha = 1;
% mu = 1;
% dirichlet = [1 1];

BETA = [1 2 4 6 10 14 20 26 34 42];
zeta = 100 * [1 1];

N = 2;

beta = BETA(N+1);

prec = 10;

err = [];

for i=2:10
    fprintf("iterate: %i\n",i);
%     simulation = [center - 2^i * l, center + 2^i * l] ;
    % mesh generation
    [Edge,E2edge,E2size,E2E,normal,K] = mesh_generation(2^i * prec,real,simulation,"mid");
    % basis function generation
    [leg_b,leg_d,dx] = basis_function(N,E2edge,Edge,E2size);
    
    % solver
    U = solver_DG_dirac_dirichlet(K,N,alpha,mu,beta,dirichlet_s,leg_b,dx,leg_d,Edge,E2edge,E2size,f,real,dirichlet,zeta);
    
    % err
    err = [err;total_error(N,real,Edge,leg_b,U,10000,sol)];
    
    figure(1)
    clf;
    view_simulation(K,N,Edge,E2edge,leg_b,U,21);
    
    figure(2)
    clf;
    view_real(K,N,real,Edge,E2edge,leg_b,U,21);
    hold on;
    fplot(sol,real,'r')
    
    drawnow;
end




order = log(err(1:end-1,:)./err(2:end,:))/log(2);
