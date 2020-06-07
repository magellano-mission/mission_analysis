function [value, isterminal, direction] = EThypPericenterEvent(t, Y, data)
%{ 
EVENT Function to stop the integration at the pericenter hyperbola;
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

% Computation of the true anomaly
vr = dot(V, R)/r;
th = acos( max( min( dot(E, R)/(e*r) , 1) , -1) );

if vr < 0
    th = - th;
end

value(1) = th;
value(2) = r - data.R_SOI;

if t > 1
    isterminal(1:2) = 1;
else 
    isterminal(1:2) = 0;
end
direction(1:2) = 0;
