function mass = mass_matrix_interface_int(K,N,real,leg_b,Edge,E2edge,E2bound,gauss)

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
        points = (x2-x1)/2 * gauss(1,:) + (x2+x1)/2;
        for i=0:N
            for j = 0:N
                f = @(x) leg_b(k,i,x).*leg_b(k,j,x);
                mass(E + i+1,E +j+1) = (x2-x1)/2 * gauss_legendre_1d(f,points,gauss(2,:));
            end
        end
    end
end

end