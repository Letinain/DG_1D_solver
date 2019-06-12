function [flux,bound] = flux_matrix_Neuman(K,N,beta,neumann,leg_b,dx,leg_d,Edge,E2edge)

ul = neumann(1);
ur = neumann(2);

bound = spalloc((N+1)*K,1,2*(N+1));
flux = spalloc((N+1)*K,(N+1)*K,3*(N+1)^2*K);

for k=1:K
    E = (k-1)*(N+1);
    
    %%% left bound
    En = (k-2)*(N+1);
    x = Edge(E2edge(k,1));
    if (k==1)
        %%%%% Bound
        for i=0:N
            bound(E + i+1) = bound(E + i+1) + 1/2 * neumann(1) * leg_b(k,i,x);
        end
        %%%%% Exterior
        for i = 0:N
            for j = 0:N
                flux(i+1,j+1) = flux(i+1,j+1) - 1/2 * leg_d(k,j,x) * leg_b(k,i,x);
            end
        end
    else
        %%%%% Interior
        for i=0:N
            for j=0:N
                flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + ( - beta / dx(k,k-1) * leg_b(k,j,x) - leg_d(k,j,x) / 2) * leg_b(k,i,x);
            end
        end
        %%%%% Exterior
        for i=0:N
            for j=0:N
                flux(E +i+1,En +j+1) = flux(E +i+1,En +j+1) + ( beta / dx(k,k-1) * leg_b(k-1,j,x) - leg_d(k-1,j,x)/2)* leg_b(k,i,x);
            end
        end
    end
    
    %%% right bound
    En = k*(N+1);
    x = Edge(E2edge(k,2));
    
    if (k==K)
        %%%%% Bound
        for i=0:N
            bound(E + i+1) = bound(E + i+1) - 1/2 * neumann(1) * leg_b(k,i,x);
        end
        %%%%% Interior
        for i = 0:N
            for j = 0:N
                flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + 1/2 * leg_d(k,j,x) * leg_b(k,i,x);
            end
        end
    else
        %%%%% Interior
        for i=0:N
            for j=0:N
                flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + ( - beta / dx(k,k+1) * leg_b(k,j,x) + leg_d(k,j,x) / 2) * leg_b(k,i,x);
            end
        end
        %%%%% Exterior
        for i=0:N
            for j=0:N
                flux(E +i+1,En +j+1) = flux(E +i+1,En +j+1) + ( beta / dx(k,k+1) * leg_b(k+1,j,x) + leg_d(k+1,j,x)/2)* leg_b(k,i,x);
            end
        end
    end
end

end