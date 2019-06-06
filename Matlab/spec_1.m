function Edge = spec_1(K,real,simulation)

a = simulation(1);
b = simulation(2);

c = real(1);
d = real(2);

l = abs(d-c);
h = l/K;

Edge_mid = [c - 0.9*l];
while (Edge_mid(end)<(d+0.9*l))
    Edge_mid = [Edge_mid,Edge_mid(end)+h];
end

Edge_left = Edge_mid(1) - 1.2*h;
i = 2;
while (Edge_left(1)>a)
    Edge_left = [Edge_left(1) - 1.2^i*h,Edge_left];
    i = i+1;
end
if (Edge_left(1)<a)
    Edge_left = Edge_left(2:end);
end
Edge_left(1) = a;

Edge_right = Edge_mid(end) + 1.2*h;
i = 2;
while (Edge_right(end)<b)
    Edge_right = [Edge_right,Edge_right(end) + 1.2^i*h];
    i = i+1;
end
if (Edge_right(end)>b)
    Edge_right = Edge_right(1:end-1);
end
if (~isempty(Edge_right)) 
    Edge_right(end) = b;
end


Edge = [Edge_left,Edge_mid,Edge_right];

end