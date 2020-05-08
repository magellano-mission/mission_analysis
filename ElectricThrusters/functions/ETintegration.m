function dY = ETintegration(t, Y, data)


M = data.M;                         % structural mass

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
    T = data.T;             % thrust [N]
else
    T = 0;
end


% Rotation Matrix
t_hat = V/norm(V);
h_hat = cross(R, V)/norm(cross(R, V));
n_hat = cross(h_hat, t_hat);
A = [t_hat, n_hat, h_hat];

aT = [T/M*1e-3; 0; 0];                          % [km/s^2]
aT_tnh = A*aT;                                  % rotated thrust

% Derivative of the states
dY(1:3) = V;
dY(4:6) = - data.mi/norm(R)^3.*R + aT_tnh;

dY = dY';

end