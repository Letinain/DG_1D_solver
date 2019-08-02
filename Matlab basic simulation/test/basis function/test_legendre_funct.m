N = 10;

LEG = zeros(N+1);
for i = 0:N
    LEG(i+1,end-i:end) = LegendrePoly(i);
end

leg_f = @(n,x) polyval(LEG(n+1,:),x);

figure(1);
clf;
hold on;
for i=0:N
    fplot(@(x) leg_f(i,x),[-1 1]);
end

test = zeros(N+1,1);
for i=0:10
    test(i+1) = integral(@(x) abs(leg_f(i,x)-legendreP(i,x)),-1,1);
end

test_logical = test <= eps;