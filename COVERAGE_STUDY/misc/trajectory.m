function [T,Y] = trajectory(X0, tspan, mu, opt)

[T, X] = ode113(@gauss, tspan, X0, opt);

n = length(T);
Y = zeros(n, 3);

for kk = 1 : n
[rr, ~] = kep2car(X(kk,:), mu);
Y(kk,:) = rr';
end

end