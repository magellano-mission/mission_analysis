function dY = EThypIntegration(t, Y, data)
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

DayActual = data.InitDay + t/86400;

KepS = uplanet(DayActual, 4);
rM2S = - kep2car_r_only(KepS);                  % retriving the sun position vector wrt the planet
light = los(R, rM2S, data.R_pl);                % checking the line of sight

if light
    T = data.T;                                 % thrust [N]
else
    T = 0;
end

T_car = OptimThrust(R, V, T, data);           

% Derivative of the states
dY(1:3) = V;
dY(4:6) = - data.mi/norm(R)^3.*R + T_car;

dY = dY';

end