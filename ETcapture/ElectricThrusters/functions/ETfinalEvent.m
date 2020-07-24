function [value, isterminal, direction] = ETfinalEvent(t, Y, data)
%{ 
EVENT: blocca quando stai sulla circolare o r = r finale dell'orbita
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
ec = norm(E);

value = ec - 2.9e-3;

if t > 1
    isterminal = 1;
else 
    isterminal = 0;
end

direction = 0;





