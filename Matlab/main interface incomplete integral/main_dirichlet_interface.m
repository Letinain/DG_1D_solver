addpath(genpath(pwd+"/.."));

% simulation domain
a = 0; b = 1;
simulation = [a b];


% real domain
c = 0; d = 0.7;
real = [c d];
% real = simulation;

[f,sol,alpha,mu,dirichlet] = sin_poisson(5,real);

BETA = [1 2 4 6 10 14 20 26 34 42];

N = 1;

beta = BETA(N+1);
err = [];

n1 = 0;
n2 = 15;

mem_edge = zeros(n2-n1+1,6);
indice = 0;

for i=n1:n2
    indice = indice +1;
    
    fprintf("iterate: %i\n",indice);
    
    % mesh generation
    [Edge,E2edge,E2size,E2E,E2bound,normal,K] = mesh_generation_interface(2^i,simulation,real,"regular");
    
    % basis function generation
    [leg_b,leg_d,dx] = basis_function(N,E2edge,Edge,E2size);
    
    % solver
    [U,A,F,mass,stiff,flux,source,bound] = solver_DG_dirichlet_interface(K,N,real,alpha,mu,beta,f,dirichlet,leg_b,leg_d,dx,Edge,E2edge,E2E,E2size,E2bound,normal);
    
    % err
    err = [err;total_error(N,real,Edge,leg_b,U,10000,sol)];
    
    mem_edge(indice,:) = [Edge(1),c,Edge(2),Edge(end-1),d,Edge(end)];
    
    figure(1);
    clf;
    view_real(K,N,real,Edge,E2edge,leg_b,U,21);
    hold on;
    fplot(sol,real,'r')
    
end

order = log(err(1:end-1,:)./err(2:end,:))/log(2);
prop_elem_bound = [(mem_edge(:,3)-mem_edge(:,2))./(mem_edge(:,3)-mem_edge(:,1))];

figure(2);
subplot(1,3,1);
plot(mem_edge(:,1:3),n1:n2);
subplot(1,3,2);
plot(order,n1+0.5:n2-0.5);
subplot(1,3,3);
plot(mem_edge(:,4:6),n1:n2);

figure(3);
loglog(err);

figure(4)
plot(log(err)/log(2));
o = zeros(1,3);
for i=1:3
    tmp = polyfit(n1:n2,-log(err(:,i)')/log(2),1);
    o(i) = tmp(1);
end