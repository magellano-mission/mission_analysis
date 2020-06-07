function dY = ETcaptIntegration(t, Y, data)
%{
ODE function to integrate the Thrusted electric arches
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

aT_car = ETcaptThrust(t, R, V, data);

% Derivative of the states
dY(1:3) = V;
dY(4:6) = - data.mi/norm(R)^3.*R + aT_car;

dY = dY';

end