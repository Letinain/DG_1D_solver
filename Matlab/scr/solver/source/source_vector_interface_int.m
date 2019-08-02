function source = source_vector_interface_int(K,N,f,leg_b,real,Edge,E2edge,E2bound,gauss)

source = zeros((N+1)*K,1);

for k=1:K
    E = (k-1)*(N+1)+1;
    if ((E2bound(k,1)==0)||(E2bound(k,2)==0))
        for i =0:N
            s = @(x) f(x).*leg_b(k,i,x);
            x1 = Edge(E2edge(k,1));
            x2 = Edge(E2edge(k,2));
            points = (x2-x1)/2 * gauss(1,:) + (x2+x1)/2;
            source(E +i) = (x2-x1)/2 * gauss_legendre_1d(s,points,gauss(2,:));
        end
    else
        if(E2bound(k,1))
            x1 = real(1);
            x2 = Edge(E2edge(k,2));
        else
            x1 = Edge(E2edge(k,1));
            x2 = real(2);
        end
        points = (x2-x1)/2 * gauss(1,:) + (x2+x1)/2;
        for i=0:N
            g = @(x) f(x).*leg_b(k,i,x);
            source(E + i) = (x2-x1)/2 * gauss_legendre_1d(g,points,gauss(2,:));
        end
    end
end

end