function [Edge,E2edge,E2size,E2E,normal,K] = mesh_division_interface(div,mesh,real)

i1 = real(1);
i2 = real(2);

K = div*(length(mesh)-1);

Edge = zeros(1,K+1);
Edge(1) = mesh(1);
for i = 2:length(mesh)
    Edge(2*(i-1)) = (mesh(i-1)+mesh(i))/2;
    Edge(2*i-1) = mesh(i);
end

ed1 = 1; ed2 = length(Edge);
while (i1 >= Edge(ed1+1))
    ed1 = ed1 + 1;
end
while (i2 <= Edge(ed2-1))
    ed2 = ed2 -1;
end

Edge = Edge(ed1:ed2);
K = length(Edge)-1;

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