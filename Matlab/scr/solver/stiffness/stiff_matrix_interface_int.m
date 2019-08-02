function stiff = stiff_matrix_interface_int(K,N,real,leg_d,Edge,E2edge,E2size,E2bound,gauss)

stiff_ref = stiff_matrix_ref(N);

stiff = spalloc((N+1)*K,(N+1)*K,(N+1)^2*K);

for k=1:K
    if((E2bound(k,1)==0)&&(E2bound(k,2)==0))
        stiff( (k-1)*(N+1)+1 : k*(N+1) , (k-1)*(N+1)+1 : k*(N+1) ) = 2 / E2size(k)^2 * stiff_ref;
    else
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
            for j=0:N
                f = @(x) leg_d(k,i,x) .* leg_d(k,j,x);
                stiff(E +i+1,E +j+1) = (x2-x1)/2 * gauss_legendre_1d(f,points,gauss(2,:));
             %   stiff(E +j+1,E +i+1) = stiff(E +i+1,E +j+1);
            end
        end
    end
end

end