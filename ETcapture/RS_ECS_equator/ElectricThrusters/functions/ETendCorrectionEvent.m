function [value, isterminal, direction] = ETendCorrectionEvent(t, Y, data)
%{ 
EVENT: blocca alla fine della manovra
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

kep = car2kep(R, V, data.mi);
kep02 = data.kep02;

value(1) = sin(kep(5) + kep(6));


if t > 1
    isterminal(1) = 1;
else 
    isterminal(1) = 0;
end

direction(1) = 0;