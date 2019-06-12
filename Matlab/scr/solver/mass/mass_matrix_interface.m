function mass = mass_matrix_interface(K,N,real,leg_b,Edge,E2edge,E2bound)

mass = speye(K*(N+1));

for k=1:K
    if ((E2bound(k,1)~=0)||(E2bound(k,2)~=0))
        E = (k-1)*(N+1);
        if(E2bound(k,1)==1)
            x1 = real(1);
            x2 = Edge(E2edge(k,2));
        else
            x1 = Edge(E2edge(k,1));
            x2 = real(2);
        end
        for i=0:N
            f = @(x) leg_b(k,i,x).^2;
            mass(E + i+1,E +i+1) = integral(f,x1,x2);
        end
    end
end

end