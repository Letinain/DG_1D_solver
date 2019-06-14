function [leg_b_n,leg_d_n,dx_n] = basis_function_interface_normalized(N,E2edge,Edge,E2size,E2bound,real)

[leg_b,leg_d,~] = basis_function(N,E2edge,Edge,E2size);

K = length(Edge)-1;

normalization = zeros(N+1,K);

for k = 1:K
    if ((E2bound(k,1)==0)&&(E2bound(k,2)==0))
        for n = 0:N
            normalization(n+1,k) = 1;
        end
    else
        if (E2bound(k,1))
            x1 = real(1);
            x2 = Edge(E2edge(k,2));
        else
            x2 = real(2);
            x1 = Edge(E2edge(k,1));
        end
        for n = 0:N
            normalization(n+1,k) = sqrt(integral(@(x) leg_b(k,n,x).^2,x1,x2));
        end
    end
end

leg_b_n = @(k,n,x) leg_b(k,n,x)/normalization(n+1,k);
leg_d_n = @(k,n,x) leg_d(k,n,x)/normalization(n+1,k);

len = Edge(2:end) - Edge(1:end-1);
len(1) = Edge(2)-real(1);
len(end) = real(2)-Edge(end-1);

dx_n = @(k1,k2) (len(k1)+len(k2))/2;

end
