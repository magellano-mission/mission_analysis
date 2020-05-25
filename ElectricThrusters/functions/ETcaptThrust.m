function aT_car = ETcaptThrust(t, R, V, M, data)

DayActual = data.InitDay + t/86400; 

KepS = uplanet(DayActual, 4);
rM2S = - kep2car_r_only(KepS);                  % retriving the sun position vector wrt the planet
light = los(R, rM2S, data.R_pl);                % checking the line of sight

if light
    T = data.T;             % thrust [N]
else
    T = 0;
end

th_range = ThRange(R, V, data);

if not(th_range)
    T = 0;
end

% Rotation Matrix tnh
t_hat = V/norm(V);
h_hat = cross(R, V)/norm(cross(R, V));
n_hat = cross(h_hat, t_hat);                % non andrebbe normalizzato? 
A = [t_hat, n_hat, h_hat];

aT_tnh = [T/M*1e-3; 0; 0];                          % [km/s^2]

aT_car = A*aT_tnh;                                  % from tnh to car

end