%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is the main script of a 1D DG
% method, with Legendre polynomials and Dirichlet
% conditions.
%
% The equation solved is:
% - alpha u'' + mu u = f
% on the interval [a,b]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
% domain
a = 0;
b = 1;


% Problem parameters
alpha = 1/2;
%mu = 1/2 * ((b-a)/pi)^2;
mu = 2/pi^2;

beta = 1;


% Dirichlet
ul = 1;
ur = 0;

% source
%f = @(x) sin(pi / (b-a) * (x+a));
f = @(x) cos(pi/2 * (x+a) / (b-a));

% mesh parameters
K = 10; % number of element
N = 30; % element degree

% subdivion of the interval
%Edge = linspace(a,b,K+1); % regular subdivision
Edge = sort([a , a + (b-a)*rand(1,K-1) , b]); % random subdivision

% element to edge matching array
E2edge = zeros(K,2);
for i = 1:K
    E2edge(i,:) = [i , i+1];
end
% Storage method [left right]

% element to edge matching array
E2E = zeros(K,2);
for i = 1:K
    E2E(i,:) = [i-1 , i+1];
end
% Storage method [left right]


% normal vector
normal = [-1,1]; % [left normal, right normal]

% elemet size array
E2size = zeros(1,K);
for i = 1:K
    E2size(i) = Edge(E2edge(i,2)) - Edge(E2edge(i,1));
end

LEG = legendre_coef(N);
LEG_d = poly_deriv(LEG);

legendre_func = @(n,x) poly_eval(n,LEG(n+1,:),x);
legendre_deriv = @(n,x) poly_eval(n-1,LEG_d(n+1,:),x);

leg_b = @(k,n,x) sqrt((2*n+1)/abs(Edge(k+1)-Edge(k))) .* legendre_func(n, 2 / abs(Edge(k+1)-Edge(k)) .* (x - (Edge(k)+Edge(k+1))/2));
leg_d = @(k,n,x) 2 / abs(Edge(k+1)-Edge(k)) .* sqrt((2*n+1)/abs(Edge(k+1)-Edge(k))) .* legendre_deriv(n, 2 / abs(Edge(k+1)-Edge(k)) * (x - (Edge(k)+Edge(k+1))/2));

%leg_b = @(k,i,x) legendre_basis(Edge(k),Edge(k+1),i,x);
%leg_d = @(k,i,x) legendre_derivative(Edge(k),Edge(k+1),i,x);
dx = @(k1,k2) (E2size(k1)+E2size(k2))/2;

%%

%%%%%%%%%%%%%
% mass matrix
%%%%%%%%%%%%%

% reference matrix
% mass_ref = eye(N+1);

% mass matrix
mass = eye((N+1)*K);
% for i = 1:K
%     mass( (i-1)*(N+1)+1 : i*(N+1) , (i-1)*(N+1)+1 : i*(N+1) ) = E2size(i)/2 * mass_ref;
% end

%%%%%%%%%%%%%%%%%%
% stiffness matrix
%%%%%%%%%%%%%%%%%%

% reference matrix
stiff_ref = zeros(N+1);

for i=0:N
    for j=0:N
        if (mod(i-j,2)==1)
            stiff_ref(i+1,j+1) = 0;
        elseif (mod(i,2)==0)
            a = min(i,j);
            stiff_ref(i+1,j+1) = sqrt((2*i+1)*(2*j+1)) * a * (a+1) ;
        else
            a = min(i,j);
            stiff_ref(i+1,j+1) = sqrt((2*i+1)*(2*j+1)) * ((a-1)*(a+2)+2);
        end
    end
end


% stiffness matrix
stiff = zeros((N+1)*K);
for i=1:K
    stiff( (i-1)*(N+1)+1 : i*(N+1) , (i-1)*(N+1)+1 : i*(N+1) ) = 2 / E2size(i)^2 * stiff_ref;
end

%%

%%%%%%%%%%%%%
% flux matrix
%%%%%%%%%%%%%

bound = zeros((N+1)*K,1);
flux = zeros((N+1)*K);

%leg_b = @(k,i,x) legendre_basis(Edge(k),Edge(k+1),i,x);
%leg_d = @(k,i,x) legendre_derivative(Edge(k),Edge(k+1),i,x);
%dx = @(k1,k2) (E2size(k1)+E2size(k2))/2;

% first cell
k = 1;
E = (k-1)*(N+1);
%%% left bound
x = Edge(1);
for i=0:N
    bound(i+1) = - beta / dx(k,k) * ul * leg_b(k,i,x);
end
for i = 0:N
    for j = 0:N
        flux(i+1,j+1) = flux(i+1,j+1) + (-beta/dx(k,k) * leg_b(k,j,x) - leg_d(k,j,x) )* leg_b(k,i,x);
    end
end
%%% right bound
En = k*(N+1);
x = Edge(E2edge(k,2));
%%%%% Interior
for i=0:N
    for j=0:N
        flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + ( - beta / dx(k,k+1) * leg_b(k,j,x) + leg_d(k,j,x) / 2) * leg_b(k,i,x);
    end
end
%%%%% Exterior
for i=0:N
    for j=0:N
        flux(E +i+1,En +j+1) = flux(E +i+1,En +j+1) + ( beta / dx(k,k+1) * leg_b(k+1,j,x) + leg_d(k+1,j,x)/2)* leg_b(k,i,x);
    end
end

% middle cell

for k=2:K-1
    E = (k-1)*(N+1);
    %%% left bound
    En = (k-2)*(N+1);
    x = Edge(E2edge(k,1));
    %%%%% Interior
    for i=0:N
        for j=0:N
            flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + ( - beta / dx(k,k-1) * leg_b(k,j,x) - leg_d(k,j,x) / 2) * leg_b(k,i,x);
        end
    end
    %%%%% Exterior
    for i=0:N
        for j=0:N
            flux(E +i+1,En +j+1) = flux(E +i+1,En +j+1) + ( beta / dx(k,k-1) * leg_b(k-1,j,x) - leg_d(k-1,j,x)/2)* leg_b(k,i,x);
        end
    end
    
    %%% right bound
    En = k*(N+1);
    x = Edge(E2edge(k,2));
    %%%%% Interior
    for i=0:N
        for j=0:N
            flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + ( - beta / dx(k,k+1) * leg_b(k,j,x) + leg_d(k,j,x) / 2) * leg_b(k,i,x);
        end
    end
    %%%%% Exterior
    for i=0:N
        for j=0:N
            flux(E +i+1,En +j+1) = flux(E +i+1,En +j+1) + ( beta / dx(k,k+1) * leg_b(k+1,j,x) + leg_d(k+1,j,x)/2)* leg_b(k,i,x);
        end
    end
end

% last cell
k = K;
E = (k-1)*(N+1);
%%% left bound
En = (k-2)*(N+1);
x = Edge(E2edge(k,1));
%%%%% Interior
for i=0:N
    for j=0:N
        flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + ( - beta / dx(k,k-1) * leg_b(k,j,x) - leg_d(k,j,x) / 2) * leg_b(k,i,x);
    end
end
%%%%% Exterior
for i=0:N
    for j=0:N
        flux(E +i+1,En +j+1) = flux(E +i+1,En +j+1) + ( beta / dx(k,k-1) * leg_b(k-1,j,x) - leg_d(k-1,j,x)/2)* leg_b(k,i,x);
    end
end
%%% right bound
x = Edge(E2edge(k,2));
for i=0:N
    bound(E + i+1) = -beta/dx(k,k) * ur * leg_b(k,i,x);
end
for i = 0:N
    for j = 0:N
        flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + (-beta/dx(k,k) * leg_b(k,j,x) + leg_d(k,j,x)) * leg_b(k,i,x);
    end
end

%%

%%%%%%%%%%%%%
% source term
%%%%%%%%%%%%%

source = zeros((N+1)*K,1);
for k = 1:K
    E = (k-1)*(N+1)+1;
    for i =0:N
        s = @(x) f(x).*leg_b(k,i,x);
        x1 = Edge(E2edge(k,1));
        x2 = Edge(E2edge(k,2));
        source(E +i) = integral(s,x1,x2);
    end
end


%%%%%%%%%%
% Assembly
%%%%%%%%%%

A = alpha * mass + mu * (stiff - flux);
F = source - mu* bound;

%%%%%%%%%%%%
% resolution
%%%%%%%%%%%%

U = A\F;


%%

%%%%%%
% view
%%%%%%
X_ref = -1:0.1:1;
l = length(X_ref);
weight = [];
for i=0:N
    weight = [weight ; sqrt((2*i+1)/2) * legendreP(i,X_ref)];
end

X = [];
Y = [];
for k = 1:K
    X = [X , E2size(k)/2 * X_ref + (Edge(E2edge(k,1))+Edge(E2edge(k,2)))/2];
    Y = [Y , sqrt(2/E2size(k)) * U((k-1)*(N+1)+1:k*(N+1))' *  weight];
end

clf;
hold on;
plot(X,Y);
%sol = @(x) sin(pi*x);
sol = @(x) cos(pi/2*x);
fplot(sol,[Edge(1) Edge(end)]);

%%%%%%%
% error
%%%%%%%
err = zeros(1,3);
err(1) = trapz(X,Y-sol(X));
err(2) = sqrt(trapz(X,(Y-sol(X)).^2));
err(3) = max(Y-sol(X));

