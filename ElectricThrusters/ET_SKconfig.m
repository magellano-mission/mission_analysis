%%
InitData = [2025 01 01 0 0 0 ];                         % initial date
data.InitDay = date2mjd2000(InitData);                  % initial day
data.TangThrust = true;                                 % true to change sma, false to change i

data.mi = astroConstants(14);
V = 1.6318e11;                                          % mars volume [km^3]
data.R_pl = nthroot(3*V/(4*pi), 3);                     % mars equivalent radius [km]

data.M = 100;                                           % [kg] Sc mass
data.sma = 6400;                                        % [km] sma
data.sma_change = 0.6;                                  % [km] sma change needed
data.i = 0*pi/180;                                      % [rad] inclination

Kep0 = [data.sma, 0, data.i , 0, 0, 0];                 % Sc initial condition
[R0, V0] = kep2car(Kep0, data.mi);
data.Y0 = [R0; V0];

data.opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Event', @ETsmaEvent);

clearvars -except data Kep0