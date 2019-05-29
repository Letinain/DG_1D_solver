function source = source_vector(K,N,f,leg_b,Edge,E2edge)

source = zeros((N+1)*K,1);
for k = 1:K
    E = (k-1)*(N+1)+1;
    for i =0:N
        s = @(x) f(x).*leg_b(k,i,x);
        x1 = Edge(E2edge(k,1));
        x2 = Edge(E2edge(k,2));
        source(E +i) = integral(s,x1,x2);
    end
end

end