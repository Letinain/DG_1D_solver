function err_by_cell = cell_error(K,N,domain,Edge,E2edge,leg_b,U,pt,sol)

err_by_cell = zeros(K,3);
for k=1:K
    X = linspace(Edge(E2edge(k,1)),Edge(E2edge(k,2)),pt);
    Y = modal_function(leg_b,N,Edge,U,X);
    err_by_cell(k,1) = trapz(X,abs(Y-sol(X)));
    err_by_cell(k,2) = sqrt(trapz(X,(Y-sol(X)).^2));
    err_by_cell(k,3) = max(abs(Y-sol(X)));
end

figure;
hold on;
for k=1:K
    plot([Edge(E2edge(k,1)) Edge(E2edge(k,2))],[err_by_cell(k,1) err_by_cell(k,1)],'b',...
        [Edge(E2edge(k,1)) Edge(E2edge(k,2))],[err_by_cell(k,2) err_by_cell(k,2)],'r',...
        [Edge(E2edge(k,1)) Edge(E2edge(k,2))],[err_by_cell(k,3) err_by_cell(k,3)],'g');
end
title("error on all the simulation");


c = domain(1);
d = domain(2);

kc = 1;
while ~((Edge(E2edge(kc,1))<= c) && (Edge(E2edge(kc,2))>= c))
    kc = kc+1;
end

kd = 1;
while ~((Edge(E2edge(kd,1))<= d) && (Edge(E2edge(kd,2))>= d))
    kd = kd+1;
end
figure;
hold on;
plot([c Edge(E2edge(kc,2))],[err_by_cell(kc,1) err_by_cell(kc,1)],'b',...
    [c Edge(E2edge(kc,2))],[err_by_cell(kc,2) err_by_cell(kc,2)],'r',...
    [c Edge(E2edge(kc,2))],[err_by_cell(kc,3) err_by_cell(kc,3)],'g');
for k=kc+1:kd-1
    plot([Edge(E2edge(k,1)) Edge(E2edge(k,2))],[err_by_cell(k,1) err_by_cell(k,1)],'b',...
        [Edge(E2edge(k,1)) Edge(E2edge(k,2))],[err_by_cell(k,2) err_by_cell(k,2)],'r',...
        [Edge(E2edge(k,1)) Edge(E2edge(k,2))],[err_by_cell(k,3) err_by_cell(k,3)],'g');
end
plot([Edge(E2edge(kd,1)) d],[err_by_cell(kd,1) err_by_cell(kd,1)],'b',...
    [Edge(E2edge(kd,1)) d],[err_by_cell(kd,2) err_by_cell(kd,2)],'r',...
    [Edge(E2edge(kd,1)) d],[err_by_cell(kd,3) err_by_cell(kd,3)],'g');
title("error on the real domain");

end
