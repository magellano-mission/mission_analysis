function [aT_car, thrust] = ETConstThrust(rM2S, R, V, M, data)

% Rotation Matrix tnh
t_hat = V/norm(V);
h_hat = cross(R, V)/norm(cross(R, V));
n_hat = cross(h_hat, t_hat)/norm(cross(h_hat, t_hat));             
A = [t_hat, n_hat, h_hat];

    d = norm(rM2S);
    [thrust] = - round(d2T(d,data), 4);
    aT_tnh = [thrust; 0; 0]/M*1e-3;             % [km/s^2]

aT_car = A*aT_tnh;                          % from tnh to car

end