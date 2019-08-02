N = 10;
E2edge = [1 2];
Edge = [-1 1];
E2size = 2;

[leg_b,leg_d,dx] = basis_function(N,E2edge,Edge,E2size);

figure(1)
hold on;
for i = 0:10
    fplot(@(x) leg_b(1,i,x),Edge);
end

figure(2)
hold on;

for i = 0:10
    fplot(@(x) leg_d(1,i,x),Edge);
end

test_integral_1 = zeros(N+1,1);
for i = 0:N
    test_integral_1(i+1) = integral(@(x) leg_b(1,i,x).^2,-1,1,'RelTol',eps,') - 1;
end

% test_integral_2 = zeros(N+1,1);
% for i = 0:N
%     test_
% end