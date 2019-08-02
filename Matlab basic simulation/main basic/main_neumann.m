% real domain
c = 0; d = 1;
real = [c d];

% n = 1;
% [f,sol,alpha,mu,dirichlet] = random_poly_poisson(n,real);
[f,sol,alpha,mu,neumann] = cos_neumann(5,real);
% n = 1;
% sol = @(x) x.^n;
% f = @(x) sol(x)-n*(n-1)*x.^(n-2);
% alpha = 1;
% mu = 1;
% dirichlet = [1 1];

BETA = [1 2 4 6 10 14 20 26 34 42];

N = 9;

beta = BETA(N+1);

err = [];

for i=0:10
    fprintf("iterate: %i\n",i);
    
    % mesh generation
    [Edge,E2edge,E2size,E2E,normal,K] = mesh_generation_basic(2^i,real,"regular");
    
    % basis function generation
    [leg_b,leg_d,dx] = basis_function(N,E2edge,Edge,E2size);
    
    % solver
    U = solver_DG_neumann(K,N,alpha,mu,beta,neumann,leg_b,dx,leg_d,Edge,E2edge,E2size,f);
    
    % err
    err = [err;total_error(N,real,Edge,leg_b,U,10000,sol)];
    clf;
    fplot(sol,real,'r');
    hold on;
    view_simulation(K,N,Edge,E2edge,leg_b,U,21);
    pause;
end

order = log(err(1:end-1,:)./err(2:end,:))/log(2);
