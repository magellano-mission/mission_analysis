function [value, isterminal, direction] = ETParab2EllipseEvent(t, Y, data)
%{ 
EVENT qui diventa parabola;
The true anomaly is computed here below and the ode stops when theta become 0.
%}

% States
R = Y(1:3);
V = Y(4:6);

% input checks
if isrow(R)
    R = R';
end

if isrow(V)
    V = V';
end

r = norm(R);
v = norm(V);
mi = data.mi;

% Computation of the eccentricity vector
E = ((v^2 - mi/r)*R - (dot(R, V)*V))/mi;
e = norm(E);

value(1) = e - 0.6;

if t > 1
    isterminal(1) = 1;
else 
    isterminal(1) = 0;
end

direction(1) = 0;