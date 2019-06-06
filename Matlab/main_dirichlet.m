% real domain
c = 0; d = 1;
real = [c d];

% n = 1;
% [f,sol,alpha,mu,dirichlet] = random_poly_poisson(n,real);
[f,sol,alpha,mu,dirichlet] = sin_poisson(real);
% n = 1;
% sol = @(x) x.^n;
% f = @(x) sol(x)-n*(n-1)*x.^(n-2);
% alpha = 1;
% mu = 1;
% dirichlet = [1 1];

BETA = [1 2 4 6 10 14 20 26 34 42];

N = 8;

beta = BETA(N+1);

err = [];
for i=0:8
    % mesh generation
    [Edge,E2edge,E2size,E2E,normal,K] = mesh_generation_basic(2^i,real,"regular");
    
    % basis function generation
    [leg_b,leg_d,dx] = basis_function(N,E2edge,Edge,E2size);
    
    % solver
    U = solver_DG_dirichlet(K,N,alpha,mu,beta,dirichlet,leg_b,dx,leg_d,Edge,E2edge,E2size,f);
    
    % err
    err = [err;total_error(N,real,Edge,leg_b,U,10000,sol)];
    
    figure(1);
    clf;
    view_real(K,N,real,Edge,E2edge,leg_b,U,21);
    hold on;
    fplot(sol,real,'r')
    
    pause;
end

order = log(err(1:end-1,:)./err(2:end,:))/log(2);
