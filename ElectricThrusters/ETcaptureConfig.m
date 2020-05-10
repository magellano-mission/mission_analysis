%%
InitData = [2025 01 01 0 0 0 ];                             % initial date
data.InitDay = date2mjd2000(InitData);                      % initial day
data.ThrustDir = 1;                                         % 1 to change sma, 2 to change e, 3 to change i

data.mi = astroConstants(14);
V = 1.6318e11;                                              % mars volume [km^3]
data.R_pl = nthroot(3*V/(4*pi), 3);                         % mars equivalent radius [km]

data.M = 7000;                                             % [kg] Sc mass
rp = 20000;
Vinf = 0.5;
e = 1 + rp*Vinf^2/data.mi;
p = rp*(e + 1);
data.R_SOI = 570000;                                             % [km] from Curtis
data.th0 = -acos((p/data.R_SOI - 1)/e)/4;
data.e_hyp = e;

data.sma = p/(data.e_hyp^2 - 1);                            % [km] sma
data.i = 25*pi/180;                                         % [rad] inclination

Kep0 = [data.sma, data.e_hyp, data.i , 0, 0, data.th0];    % Sc initial condition
[R0, V0] = kep2car(Kep0, data.mi);
data.Y0 = [R0; V0];

data.opt1 = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);% 'Event', @EThyp2ParabEvent);
data.opt2 = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);% 'Event', @EThyp2ParabEvent);
data.opt3 = odeset('RelTol', 1e-10, 'AbsTol', 1e-10); %, 'Event', @EThyp2ParabEvent);

clearvars -except data Kep0