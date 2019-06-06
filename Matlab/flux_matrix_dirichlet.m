function [flux,bound] = flux_matrix_dirichlet(K,N,beta,dirichlet,leg_b,dx,leg_d,Edge,E2edge)

ul = dirichlet(1);
ur = dirichlet(2);

bound = zeros((N+1)*K,1);
flux = zeros((N+1)*K);

for k=1:K
    E = (k-1)*(N+1);
    
    %%% left bound
    En = (k-2)*(N+1);
    x = Edge(E2edge(k,1));
    if (k==1)
        %%%%% Bound
        for i=0:N
            bound(E + i+1) = bound(E + i+1) + beta / dx(k,k) * ul * leg_b(k,i,x);
        end
        %%%%% Exterior
        for i = 0:N
            for j = 0:N
                flux(i+1,j+1) = flux(i+1,j+1) + (-beta/dx(k,k) * leg_b(k,j,x) - leg_d(k,j,x) )* leg_b(k,i,x);
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
            bound(E + i+1) = bound(E + i+1) + beta/dx(k,k) * ur * leg_b(k,i,x);
        end
        %%%%% Interior
        for i = 0:N
            for j = 0:N
                flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + (-beta/dx(k,k) * leg_b(k,j,x) + leg_d(k,j,x)) * leg_b(k,i,x);
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