function source = source_vector_interface(K,N,f,leg_b,real,Edge,E2edge,E2bound)

source = zeros((N+1)*K,1);

for k=1:K
    E = (k-1)*(N+1)+1;
    if ((E2bound(k,1)==0)||(E2bound(k,2)==0))
        for i =0:N
            s = @(x) f(x).*leg_b(k,i,x);
            x1 = Edge(E2edge(k,1));
            x2 = Edge(E2edge(k,2));
            source(E +i) = integral(s,x1,x2);
        end
    else
        if(E2bound(k,1))
            x1 = real(1);
            x2 = Edge(E2edge(k,2));
        else
            x1 = Edge(E2edge(k,1));
            x2 = real(2);
        end
        for i=0:N
            g = @(x) f(x).*leg_b(k,i,x);
            source(E + i) = integral(g,x1,x2);
        end
    end
end
    
end