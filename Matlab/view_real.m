function view_real(~,N,domain,Edge,E2edge,leg_b,U,pt)

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

if (kc==kd)
    X = linspace(c,d,pt);
else
    X = zeros((kd-kc+1)*pt,1);
    X(1:pt) = linspace(Edge(E2edge(kc,1)),Edge(E2edge(kd,2)),pt);
    for i=2:kd-kc
        X((i-1)*pt+1:i*pt) = linspace(Edge(E2edge(kc+i-1,1)),Edge(E2edge(kc+i-1,2)),pt);
    end
    X((kd-kc)*pt+1:(kd-kc+1)*pt) = linspace(Edge(E2edge(kd,1)),d,pt);
end
Y = modal_function(leg_b,N,Edge,U,X);
figure;
plot(X,Y);

end