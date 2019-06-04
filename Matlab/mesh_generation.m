function [Edge,E2edge,E2size,E2E,normal,K] = mesh_generation(prec,real,simulation,method)

a = simulation(1);
b = simulation(2);

% subdivion of the interval
if (method == "regular")
    Edge = linspace(a,b,prec+1); % regular subdivision
elseif ( method == "random")
    Edge = sort([a , a + (b-a)*rand(1,prec-1) , b]); % random subdivision
elseif ( method == "mid")
    Edge = spec_1(prec,real,simulation);
else
    error('Bad third argument. It can be "regular" or "random"');
end

K = length(Edge)-1;

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
    E2size(i) = abs(Edge(E2edge(i,2)) - Edge(E2edge(i,1)));
end


end