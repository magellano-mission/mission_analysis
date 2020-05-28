function aT_car = ETConstThrust(R, V, M, data)

% Rotation Matrix tnh
t_hat = V/norm(V);
h_hat = cross(R, V)/norm(cross(R, V));
n_hat = cross(h_hat, t_hat)/norm(cross(h_hat, t_hat));                % non andrebbe normalizzato? 
A = [t_hat, n_hat, h_hat];

if data.direction == "tangential"
    aT_tnh = [data.T; 0; 0]/M*1e-3;             % [km/s^2]
elseif data.direction == "perpendicular"
    aT_tnh = [0; 0; data.T]/M*1e-3;
end

aT_car = A*aT_tnh;                          % from tnh to car

end