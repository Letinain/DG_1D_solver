function stiff_ref = stiff_matrix_ref(N)

stiff_ref = zeros(N+1,N+1);

for i=0:N
    for j=0:N
        if (mod(i-j,2)==1)
            stiff_ref(i+1,j+1) = 0;
        elseif (mod(i,2)==0)
            minimum = min(i,j);
            stiff_ref(i+1,j+1) = sqrt((2*i+1)*(2*j+1)) * minimum * (minimum+1) ;
        else
            minimum = min(i,j);
            stiff_ref(i+1,j+1) = sqrt((2*i+1)*(2*j+1)) * ((minimum-1)*(minimum+2)+2);
        end
    end
end

end