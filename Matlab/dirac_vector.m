function [Dirac_mat,Dirac_vec] = dirac_term(K,N,domain,dirichlet,zeta,E2edge,Edge,leg_b)

Dirac_mat = zeros((N+1)*K);
Dirac_vec = zeros((N+1)*K,1);

c = domain(1);
d = domain(2);

kc = 1;
while ~((Edge(E2edge(kc,1))<= c) && (Edge(E2edge(kc,2))>= c))
    kc = kc+1;
end

for i=0:N
    for j=0:N
        Dirac_mat((kc-1)*(N+1)+i+1,(kc-1)*(N+1)+j+1) = zeta(1) * leg_b(kc,i,c) * leg_b(kc,j,c);
    end
    Dirac_vec((kc-1)*(N+1)+i+1) = dirichlet(1) * leg_b(kc,i,c);
end

kd = 1;
while ~((Edge(E2edge(kd,1))<= d) && (Edge(E2edge(kd,2))>= d))
    kd = kd+1;
end

for i=0:N
    for j=0:N
        Dirac_mat((kd-1)*(N+1)+i+1,(kd-1)*(N+1)+j+1) = zeta(2) * leg_b(kd,i,d) * leg_b(kd,j,d);
    end
    Dirac_vec((kc-1)*(N+1)+i+1) = dirichlet(2) * leg_b(kd,i,d);
end

end