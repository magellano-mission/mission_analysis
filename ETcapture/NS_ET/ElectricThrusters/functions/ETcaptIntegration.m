function dY = ETcaptIntegration(t, Y, data)
%{
ODE function to integrate the Thrusted electric arches
%}
t
% States
R = Y(1:3);
V = Y(4:6);
M = Y(7);

% input checks
if isrow(R)
    R = R';
end
if isrow(V)
    V = V';
end

aJ2_car = ETcaptGrav(t, R, data);
aT_car = ETcaptThrust(t, R, V, M, data);

% Derivative of the states
dY(1:3) = V;
dY(4:6) = - data.mi/norm(R)^3.*R + aT_car + aJ2_car;
dY(7) = - norm(data.T)/data.Isp/data.g0;

dY = dY';

end