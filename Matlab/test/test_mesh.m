f = {};
f{1} = [0 1];
for i = 1:10
[Edge,E2edge,E2size,E2E,Edge2bound,normal,K] = mesh_division(2,f{i});
f{i+1} = Edge;
end

e = {};
for i = 0:10
[Edge,E2edge,E2size,E2E,Edge2bound,normal,K] = mesh_generation_basic(2^i,[0 1],"regular");
e{i+1} = Edge;
end

%%

f = {};
f{1} = [0 1];
for i = 1:10
[Edge,E2edge,E2size,E2E,normal,K] = mesh_division_interface(2,f{i},[0 0.7]);
f{i+1} = Edge;
end

e = {};
for i = 0:10
[Edge,E2edge,E2size,E2E,E2bound,normal,K] = mesh_generation_interface(2^i,[0 1],[0 0.7],"regular");
e{i+1} = Edge;
end


%%

f = {};
f{1} = mesh_generation_interface_bb(1,[0 1],[0 0.5000001],[1 1.2]);
for i = 1:10
[Edge,E2edge,E2size,E2E,E2bound,normal,K] = mesh_division_interface_bb(2,f{i},[0 0.7],[1 1.2]);
f{i+1} = Edge;
end

e = {};
for i = 0:10
[Edge,E2edge,E2size,E2E,E2bound,normal,K] = mesh_generation_interface_bb(2^i,[0 1],[0 0.5000000001],[1 1.2]);
e{i+1} = Edge;
end

