addpath(genpath(pwd+"/.."));

BETA = [1 2 4 6 10 14 20 26 34 42];

% simulation domain
a = 0; b = 1;
simulation = [a b];

tol = 1e-5;
pas = 0.9;

n2 = floor(log(tol)/log(pas))+1;

prec = 1024;

N = 6;
beta = BETA(N+1);

mem_edge = zeros(n2+1,3);
prop = pas.^(0:n2);
err_L1 = zeros(n2+1,5);
err_L2 = zeros(n2+1,5);
err_max = zeros(n2+1,5);
cond = zeros(n2+1,5);

% basic
parfor i = 1:n2+1
    [Edge,E2edge,E2size,E2E,E2bound,normal,K,real] = mesh_generation_interface_adapt(prec,simulation,prop(i));
    [f,sol,alpha,mu,dirichlet] = sin_poisson(5,real);
    [leg_b,leg_d,dx] = basis_function_interface(N,E2edge,Edge,E2size,real);
    [U,A,~] = solver_DG_dirichlet_interface(K,N,real,alpha,mu,beta,f,dirichlet,leg_b,leg_d,dx,Edge,E2edge,E2E,E2size,E2bound,normal);
    err = total_error(N,real,Edge,leg_b,U,10000,sol);
    
    err_L1(i,1) = err(1);
    err_L2(i,1) = err(2);
    err_max(i,1) = err(3);
    cond(i,1) = condest(A);
end

parfor i = 1:n2+1
    % normalized
    [Edge,E2edge,E2size,E2E,E2bound,normal,K,real] = mesh_generation_interface_adapt(prec,simulation,prop(i));
    [f,sol,alpha,mu,dirichlet] = sin_poisson(5,real);
    [leg_b,leg_d,dx] = basis_function_interface_normalized(N,E2edge,Edge,E2size,E2bound,real);
    [U,A,~] = solver_DG_dirichlet_interface(K,N,real,alpha,mu,beta,f,dirichlet,leg_b,leg_d,dx,Edge,E2edge,E2E,E2size,E2bound,normal);
    err = total_error(N,real,Edge,leg_b,U,10000,sol);
    
    err_L1(i,2) = err(1);
    err_L2(i,2) = err(2);
    err_max(i,2) = err(3);
    cond(i,2) = condest(A);
end

parfor i = 1:n2+1
    % bounding box 100
    [Edge,E2edge,E2size,E2E,E2bound,normal,K,real] = mesh_generation_interface_bb_adapt(prec,simulation,prop(i),1);
    [f,sol,alpha,mu,dirichlet] = sin_poisson(5,real);
    [leg_b,leg_d,dx] = basis_function_interface(N,E2edge,Edge,E2size,real);
    [U,A,~] = solver_DG_dirichlet_interface(K,N,real,alpha,mu,beta,f,dirichlet,leg_b,leg_d,dx,Edge,E2edge,E2E,E2size,E2bound,normal);
    err = total_error(N,real,Edge,leg_b,U,10000,sol);
    
    err_L1(i,3) = err(1);
    err_L2(i,3) = err(2);
    err_max(i,3) = err(3);
    cond(i,3) = condest(A);
    
end

parfor i = 1:n2+1
    
    
    % bounding box 120
    [Edge,E2edge,E2size,E2E,E2bound,normal,K,real] = mesh_generation_interface_bb_adapt(prec,simulation,prop(i),1.2);
    [f,sol,alpha,mu,dirichlet] = sin_poisson(5,real);
    [leg_b,leg_d,dx] = basis_function_interface(N,E2edge,Edge,E2size,real);
    [U,A,~] = solver_DG_dirichlet_interface(K,N,real,alpha,mu,beta,f,dirichlet,leg_b,leg_d,dx,Edge,E2edge,E2E,E2size,E2bound,normal);
    err = total_error(N,real,Edge,leg_b,U,10000,sol);
    
    err_L1(i,4) = err(1);
    err_L2(i,4) = err(2);
    err_max(i,4) = err(3);
    cond(i,4) = condest(A);
end

parfor i = 1:n2+1
    % bounding box 120 + normalized
    [Edge,E2edge,E2size,E2E,E2bound,normal,K,real] = mesh_generation_interface_bb_adapt(prec,simulation,prop(i),1.2);
    [f,sol,alpha,mu,dirichlet] = sin_poisson(5,real);
    [leg_b,leg_d,dx] = basis_function_interface_normalized(N,E2edge,Edge,E2size,E2bound,real);
    [U,A,~] = solver_DG_dirichlet_interface(K,N,real,alpha,mu,beta,f,dirichlet,leg_b,leg_d,dx,Edge,E2edge,E2E,E2size,E2bound,normal);
    err = total_error(N,real,Edge,leg_b,U,10000,sol);
    
    err_L1(i,5) = err(1);
    err_L2(i,5) = err(2);
    err_max(i,5) = err(3);
    cond(i,5) = condest(A);
end

figure;
subplot(2,2,1);
p = semilogy(prop,err_L1);
xlabel('ratio');
ylabel('L1 error');
legend('basic','normalized','bb 100','bb 120','normalized bb 120');
p(1).Marker = '+';
p(4).Marker = '*';

subplot(2,2,2);
p = semilogy(prop,err_L2);
xlabel('ratio');
ylabel('L2 error');
p(1).Marker = '+';
p(4).Marker = '*';

subplot(2,2,3);
p = semilogy(prop,err_max);
xlabel('ratio');
ylabel('max error');
p(1).Marker = '+';
p(4).Marker = '*';

subplot(2,2,4);
p = semilogy(prop,cond);
xlabel('ratio');
ylabel('matrix condition number');
p(1).Marker = '+';
p(4).Marker = '*';
