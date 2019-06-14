addpath(genpath(pwd+"/.."));

BETA = [1 2 4 6 10 14 20 26 34 42];

% simulation domain
a = 0; b = 1;
simulation = [a b];

tol = 1e-5;
pas = 0.9;

n2 = floor(log(tol)/log(pas))+1;

N = 3;
beta = BETA(N+1);

prec = 2.^5;

mem_edge = zeros(n2+1,3);
err = zeros(n2+1,3);
cond = zeros(n2+1,1);

parfor i = 0:n2
    
    [Edge,E2edge,E2size,E2E,E2bound,normal,K,real] = mesh_generation_interface_adapt(prec,simulation,pas^i);
    
    [f,sol,alpha,mu,dirichlet] = sin_poisson(5,real);
    
    [leg_b,leg_d,dx] = basis_function_interface_normalized(N,E2edge,Edge,E2size,E2bound,real);
    
    % solver
    [U,A,F] = solver_DG_dirichlet_interface(K,N,real,alpha,mu,beta,f,dirichlet,leg_b,leg_d,dx,Edge,E2edge,E2E,E2size,E2bound,normal);
    
    
    
    % err
    err(i+1,:) = total_error(N,real,Edge,leg_b,U,10000,sol);
    
    mem_edge(i+1,:) = [Edge(end-1),real(2),Edge(end)];
    
    cond(i+1) = condest(A);
    
    % figure(1);
    % clf;
    % view_real(K,N,real,Edge,E2edge,leg_b,U,21);
    % hold on;
    % fplot(sol,real,'r')
end

e1 = mem_edge(1,1);
e2 = mem_edge(1,3);

fig = figure(1);
clf;

subplot(3,2,1);
title(sprintf('N = %d, prec = %d, dx = %.2e',N,prec,1/prec));
plot(err);
xlabel("iteration");
ylabel("error");

subplot(3,2,3);
plot((mem_edge(:,2)-e1)/(e2-e1),log(cond)/log(10));
xlabel("ratio C-/C");
ylabel("matrix condition number");

subplot(3,2,2);
plot(mem_edge(:,2),err);
xlabel("interface x position");
ylabel("error");

subplot(3,2,4);
plot((mem_edge(:,2)-e1)/(e2-e1),err);
xlabel("ratio C-/C");
ylabel("error");

subplot(3,2,5);
plot(0:n2,mem_edge);
xlabel("iteration");
ylabel("interface x position");

subplot(3,2,6);
plot(0:n2,(mem_edge(:,2)-e1)/(e2-e1));
xlabel("iteration");
ylabel("ratio C-/C");


% calcul prop de mauvis conditionnement
n_max_d = floor(n2/2);
lgcond = log(cond)/log(10);
D2lgcond = lgcond(3:n_max_d) + lgcond(1:n_max_d-2) - 2*lgcond(2:n_max_d-1);
[~,ind] = max(D2lgcond);
cut_ratio_cond = (mem_edge(ind+1,2)-e1)/(e2-e1);
