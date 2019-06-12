function [flux,bound] = flux_matrix_seconde(K,N,beta,seconde,leg_b,dx,leg_d,leg_dd,Edge,E2edge)

ul = seconde(1);
ur = seconde(2);

beta_0 = beta(1);
beta_1 = beta(2);

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
            bound(E + i+1) = bound(E + i+1) - beta_1 * dx(k,k) * ul * leg_b(k,i,x);
        end
        %%%%% Exterior
        for i = 0:N
            for j = 0:N
                flux(i+1,j+1) = flux(i+1,j+1) + (- leg_d(k,j,x) - beta_1 * dx(k,k) * leg_dd(k,i,x)) * leg_b(k,i,x);
            end
        end
    else
        %%%%% Interior
        for i=0:N
            for j=0:N
                flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + ( - beta_0 / dx(k,k-1) * leg_b(k,j,x) - leg_d(k,j,x) / 2 - beta_1 * dx(k,k-1) * leg_dd(k,j,x) ) * leg_b(k,i,x);
            end
        end
        %%%%% Exterior
        for i=0:N
            for j=0:N
                flux(E +i+1,En +j+1) = flux(E +i+1,En +j+1) + ( beta_0 / dx(k,k-1) * leg_b(k-1,j,x) - leg_d(k-1,j,x)/2 + beta_1 * dx(k,k-1) * leg_dd(k-1,j,x) )* leg_b(k,i,x);
            end
        end
    end
    
    %%% right bound
    En = k*(N+1);
    x = Edge(E2edge(k,2));
    
    if (k==K)
        %%%%% Bound
        for i=0:N
            bound(E + i+1) = bound(E + i+1) - beta_1 * dx(k,k) * ur * leg_b(k,i,x);
        end
        %%%%% Interior
        for i = 0:N
            for j = 0:N
                flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1)  + (- leg_d(k,j,x) - beta_1 * dx(k,k) * leg_dd(k,i,x)) * leg_b(k,i,x);
            end
        end
    else
        %%%%% Interior
        for i=0:N
            for j=0:N
                flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + ( - beta_0 / dx(k,k+1) * leg_b(k,j,x) + leg_d(k,j,x) / 2 - beta_1 * dx(k,k+1) * leg_dd(k,j,x) ) * leg_b(k,i,x);
            end
        end
        %%%%% Exterior
        for i=0:N
            for j=0:N
                flux(E +i+1,En +j+1) = flux(E +i+1,En +j+1) + ( beta_0 / dx(k,k+1) * leg_b(k+1,j,x) + leg_d(k+1,j,x)/2 + beta_1 * dx(k,k+1) * leg_dd(k+1,j,x) )* leg_b(k,i,x);
            end
        end
    end
end

end