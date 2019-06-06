function [Edge,E2edge,E2size,E2E,normal,K] = mesh_generation_basic(prec,real,mod)

a = real(1);
b = real(2);

if (mod == "regular")
    Edge = linspace(a,b,prec+1);
elseif (mod == "random")
    Edge = [a , sort(rand(1,prec-1)) , b];
else
    error("last argument must be regular or random");
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
