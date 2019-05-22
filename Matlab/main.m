%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is the main script of a 1D DG
% method, with Legendre polynomials and Dirichlet
% conditions.
% 
% The equation solved is:
% - alpha u'' + mu u = f
% on the interval [a,b]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%global alpha mu K a b N Edge E2edge normal E2size mass;

% Problem parameters
alpha = 1;
mu = 1;

beta = 1;

a = 5;
b = 20;

K = 10; % number of element
N = 5; % element degree

% subdivion of the interval
% Edge = linspace(a,b,K+1); % regular subdivision
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

%%%%%%%%%%%%%
% mass matrix
%%%%%%%%%%%%%

% reference matrix
mass_ref = eye(N+1);

% mass matrix
mass = zeros((N+1)*K);
for i = 1:K
    mass( (i-1)*(N+1)+1 : i*(N+1) , (i-1)*(N+1)+1 : i*(N+1) ) = E2size(i)/2 * mass_ref;
end

%%%%%%%%%%%%%%%%%%
% stiffness matrix
%%%%%%%%%%%%%%%%%%

% reference matrix
stiff_ref = zeros(N+1);

for i=1:N+1
    for j=1:N+1
        if (mod(i-j,2)==1)
            stiff_ref(i,j) = 0;
        elseif (mod(i,2)==0)
            stiff_ref(i,j) = sqrt((2*i+1)*(2*j+1)) * min(i,j) * (2*min(i,j)+1) ;
        else
            stiff_ref(i,j) = sqrt((2*i+1)*(2*j+1)) * (min(i,j) * (2*min(i,j)+1) + 1);
        end
    end
end


% stiffness matrix
stiff = zeros((N+1)*K);
for i=1:K
    stiff( (i-1)*(N+1)+1 : i*(N+1) , (i-1)*(N+1)+1 : i*(N+1) ) = 4 / E2size(i) * stiff_ref;
end

%%%%%%%%%%%%%
% flux matrix
%%%%%%%%%%%%%

flux = zeros((N+1)*K);

% right side matrix

% self interation flux



