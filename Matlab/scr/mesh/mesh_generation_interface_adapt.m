function [Edge,E2edge,E2size,E2E,E2bound,normal,K,real] = mesh_generation_interface_adapt(prec,simulation,prop)

s1 = simulation(1);
s2 = simulation(2);

% edge
Edge = linspace(s1,s2,prec+1);

% real domain
tmp = Edge(end-1);
real = [s1,tmp+prop*(s2 - tmp)];

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

E2bound = zeros(K,2);
E2bound(1,1) = 1;
E2bound(end,end) = 1;


end