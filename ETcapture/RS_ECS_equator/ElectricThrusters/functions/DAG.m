function aT_car = DAG(R, V, M, data)

sma_t = data.sma_t;
e_t = data.e_t;
i_t = data.i_t;
RAAN_t = data.RAAN_t;

Y02 = data.Y02;
kep2 = car2kep(Y02(1:3), Y02(4:6), data.mi);
sma_0 = kep2(1);
e_0 = kep2(2);
i_0 = kep2(3);
RAAN_0 = kep2(4);

%% Scrivere alpha e beta per i parametri selezionati
kep = car2kep(R, V, data.mi);

% SMA
alp_sma = atan((kep(2)*sin(kep(6))/(1+kep(2)*cos(kep(6)))));
beta_sma = 0;

% Eccentricity
E = 2*atan(sqrt((1-kep(2))/(1+kep(2)))*tan(kep(6)/2));
alp_e = atan(sin(kep(6))/(cos(kep(6) + cos(E))));
beta_e = 0;

% Inclination
alp_i = 0;
beta_i = sign(cos(kep(5)+kep(6)))*(pi/2);

% RAAN
alp_RAAN = 0;
beta_RAAN = sign(cos(kep(5)+kep(6)))*(pi/2);

%% thrust directions per i parametri selezionati
f_sma = [cos(beta_sma)*cos(alp_sma); cos(beta_sma)*sin(alp_sma); sin(beta_sma)];
f_e = [cos(beta_e)*cos(alp_e); cos(beta_e)*sin(alp_e); sin(beta_e)];
f_i = [cos(beta_i)*cos(alp_i); cos(beta_i)*sin(alp_i); sin(beta_i)];
f_RAAN = [cos(beta_RAAN)*cos(alp_RAAN); cos(beta_RAAN)*sin(alp_RAAN); sin(beta_RAAN)];

%% adaptive ratio...
AR_sma = (sma_t - kep(1))/(sma_t - sma_0);
AR_e = (e_t - kep(2))/(e_t - e_0);
AR_i = (i_t - kep(3))/(i_t - i_0);
AR_RAAN = (RAAN_t - kep(4))/(RAAN_t - RAAN_0);

% Kronecker delta
if abs(kep(1) - sma_t) < 1
    delta_sma = 1;
else
    delta_sma = 0;
end

if abs(kep(2) - e_t) < 1e-3
    delta_e = 1;
else
    delta_e = 0;
end

if (kep(3) - i_t) < 1e-3
    delta_i = 1;
else
    delta_i = 0;
end

if (kep(4) - RAAN_t) < 1e-3
    delta_RAAN = 1;
else
    delta_RAAN = 0;
end

%% thrust vector nel frame tnh 
fT_vers = (1 - delta_RAAN)*AR_RAAN*f_RAAN + (1 - delta_i)*AR_i*f_i +...
     (1 - delta_sma)*AR_sma*f_sma + (1 - delta_e)*AR_e*f_e;

fT = fT_vers*data.T;     

% Rotation Matrix tnh
t_hat = V/norm(V);
h_hat = cross(R, V)/norm(cross(R, V));
n_hat = cross(h_hat, t_hat);                % non andrebbe normalizzato? 
A = [t_hat, n_hat, h_hat];

aT_tnh = fT./M*1e-3;                      % [km/s^2]

aT_car = A*aT_tnh;                          % from tnh to car

end