function stiff = stiff_matrix(K,N,E2size)

stiff_ref = stiff_matrix_ref(N);

stiff = zeros((N+1)*K);
for i=1:K
    stiff( (i-1)*(N+1)+1 : i*(N+1) , (i-1)*(N+1)+1 : i*(N+1) ) = 2 / E2size(i)^2 * stiff_ref;
end

end