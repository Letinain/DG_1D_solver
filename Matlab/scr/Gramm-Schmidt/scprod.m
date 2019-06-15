function s = scprod(fun1,fun2,x0,x1)
    s = integral(@(x) fun1(x).*fun2(x),x0,x1);
end