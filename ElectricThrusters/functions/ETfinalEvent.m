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

value(1) = ec - 5e-4;
value(2) = r - data.r_final;

if t > 1
    isterminal(1:2) = 1;
else 
    isterminal(1:2) = 0;
end

direction(1:2) = 0;





