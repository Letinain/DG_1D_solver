function Dirac = dirac_vector(K,N,domain,zeta,E2edge,Edge,leg_b)

Dirac = zeros((N+1)*K);

c = domain(1);
d = domain(2);

kc = 1;
while ~((Edge(E2edge(kc,1))<= c) && (Edge(E2edge(kc,2))>= c))
    kc = kc+1;
end

for i=0:N
    for j=0:N
        Dirac((kc-1)*(N+1)+i+1,(kc-1)*(N+1)+j+1) = zeta(1) * leg_b(kc,i,c) * leg_b(kc,j,c);
    end
end

kd = 1;
while ~((Edge(E2edge(kd,1))<= d) && (Edge(E2edge(kd,2))>= d))
    kd = kd+1;
end

for i=0:N
    for j=0:N
        Dirac((kd-1)*(N+1)+i+1,(kd-1)*(N+1)+j+1) = zeta(2) * leg_b(kd,i,d) * leg_b(kd,j,d);
    end
end

end