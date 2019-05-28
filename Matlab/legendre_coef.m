function Leg = legendre_coef(N)
Leg = zeros(N+1);

Leg(1,1) = 1;
Leg(2,2) = 1;

for k = 2:N
    Leg(k+1,:) = ((2*k-1)*[0,Leg(k,1:end-1)]-(k-1)*Leg(k-1,:))/k;
end

end