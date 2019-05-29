function view_simulation(K,N,Edge,E2edge,leg_b,U,pt)

X = zeros(K*pt,1);
for k=1:K
    X((k-1)*pt+1:k*pt) = linspace(Edge(E2edge(k,1)),Edge(E2edge(k,2)),pt);
end

Y = modal_function(leg_b,N,Edge,U,X);

figure;
plot(X,Y);
end