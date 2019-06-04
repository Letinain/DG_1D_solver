function err = total_error(N,domain,Edge,leg_b,U,pt,sol)

c = domain(1);
d = domain(2);

X = linspace(c,d,pt);
Y = modal_function(leg_b,N,Edge,U,X);
err = zeros(3,2);
err(1,1) = trapz(X,abs(Y-sol(X)));
err(1,2) = 100 * err(1,1)/trapz(X,abs(sol(X)));
err(2,1) = sqrt(trapz(X,(Y-sol(X)).^2));
err(2,2) = 100 * err(2,1)/sqrt(trapz(X,sol(X).^2));
err(3,1) = max(abs(Y-sol(X)));
err(3,2) = 100 * err(3,1)/max(abs(sol(X)));
fprintf('error L1: %e (%2.6f %%) \nerror L2: %e (%2.6f %%)\nerror max: %e (%2.6f %%)\n',err(1,1),err(1,2),err(2,1),err(2,2),err(3,1),err(3,2)); 

err = err(:,1)';

end