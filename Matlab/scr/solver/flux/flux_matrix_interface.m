function [flux,bound] = flux_matrix_interface(K,N,real,beta,dirichlet,leg_b,dx,leg_d,Edge,E2edge,E2E,E2bound,normal)

bound = spalloc((N+1)*K,1,2*(N+1));
flux = spalloc((N+1)*K,(N+1)*K,3*(N+1)^2*K);

% E2bound = @(k,n) Edge2bound(E2edge(k,n));

for k1 = 1:K
    
    E = (k1-1)*(N+1);
    
    for edge=1:2
        if(E2bound(k1,edge)==1)
            x = real(edge);
            n = normal(edge);
            %%% bound
            dir = dirichlet(edge);
            for i=0:N
                bound(E + i+1) = bound(E + i+1) + dirichlet_flux_bound(beta,leg_b,dx,k1,i,x,dir);
            end
            %%% interior
            for i = 0:N
                for j = 0:N
                    flux(E + i+1,E + j+1) = flux(E + i+1, E + j+1) + dirichlet_flux_interior(beta,leg_b,leg_d,dx,k1,i,j,x,n);
                end
            end
            
        else
            x = Edge(E2edge(k1,edge));
            n = normal(edge);
            
            k2 = E2E(k1,edge);
            En = (k2-1)*(N+1);
            
            %%%%% Interior
            for i=0:N
                for j=0:N
                    flux(E +i+1,E +j+1) = flux(E +i+1,E +j+1) + basic_flux_interior(beta,leg_b,leg_d,dx,k1,k2,i,j,x,n);
                end
            end
            
            %%%%% Exterior
            for i=0:N
                for j=0:N
                    flux(E +i+1,En +j+1) = flux(E +i+1,En +j+1) + basic_flux_exterior(beta,leg_b,leg_d,dx,k1,k2,i,j,x,n);
                end
            end
        end
        
    end
    
end

end