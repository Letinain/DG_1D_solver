function [Edge,E2edge,E2size,E2E,normal,K] = mesh_division(div,mesh)

K = div*(length(mesh)-1);

Edge = zeros(1,K+1);
Edge(1) = mesh(1);
for i = 2:length(mesh)
    Edge(2*(i-1)) = (mesh(i-1)+mesh(i))/2;
    Edge(2*i-1) = mesh(i);
end

E2edge = zeros(K,2);
for i = 1:K
    E2edge(i,:) = [i , i+1];
end

E2E = zeros(K,2);
for i = 1:K
    E2E(i,:) = [i-1 , i+1];
end

normal = [-1,1];

E2size = zeros(1,K);
for i = 1:K
    E2size(i) = abs(Edge(E2edge(i,2)) - Edge(E2edge(i,1)));
end


end