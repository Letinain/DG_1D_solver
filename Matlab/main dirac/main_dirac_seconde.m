% simulation domain
a = -1.5434361; b= 2.48414;
simulation = [a b];
seconde = [0 0];

% real domain
c = 0; d = 1;
real = [c d];

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
zeta = 1e10 * [1 1];

N = 3;

beta = [BETA(N+1) 100];

err = [];

for i=2:10
    % mesh generation
    [Edge,E2edge,E2size,E2E,normal,K] = mesh_generation(2^i,real,simulation,"regular");
    % basis function generation
    [leg_b,leg_d,leg_dd,dx] = basis_function_seconde(N,E2edge,Edge,E2size);
    
    % solver
    U = solver_DG_dirac_seconde(K,N,alpha,mu,beta,seconde,leg_b,dx,leg_d,leg_dd,Edge,E2edge,E2size,f,real,dirichlet,zeta);
    
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
    
    figure(1)
    clf;
    view_simulation(K,N,Edge,E2edge,leg_b,U,21);
    
    figure(2)
    clf;
    view_real(K,N,real,Edge,E2edge,leg_b,U,21);
    hold on;
    fplot(sol,real,'r')
    
    drawnow;
    pause;
end

order = log(err(1:end-1,:)./err(2:end,:))/log(2);
