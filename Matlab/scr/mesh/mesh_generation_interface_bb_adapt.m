function [Edge,E2edge,E2size,E2E,E2bound,normal,K,real] = mesh_generation_interface_bb_adapt(prec,simulation,prop,bb)

s1 = simulation(1);
s2 = simulation(2);

Edge = linspace(s1,s2,prec+1);

i2 = Edge(end-1) + prop * (Edge(end)-Edge(end-1));

real = [s1 i2];

Edge(end) = min(Edge(end-1) + bb*(i2-Edge(end-1)),Edge(end));

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
