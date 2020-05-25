%%
InitData = [2027 04 25 0 0 0 ];                             % initial date
data.InitDay = date2mjd2000(InitData);                      % initial day
data.ThrustDir = 1;                                         % 1 to change sma, 2 to change e, 3 to change i
data.Isp = 4300; % alternative 4300 s
data.g0 = 9.807;

data.mi = astroConstants(14);
data.switchers.grav = true;
data.J2 = 0.00196; %Mars J2
data.J3 = 3.1450e-5; %Mars J3
data.J4 = -1.5377e-5; %Mars J4
T_lla2 = 1.02749125 * 24 * 3600;   % period [s]  
data.wM = 2 * pi / T_lla2;         % [rad/s]
V = 1.6318e11;                                              % mars volume [km^3]
data.R_pl = nthroot(3*V/(4*pi), 3);                         % mars equivalent radius [km]
data.r_final = 12300; 

data.M = 2400;                                              % [kg] Sc mass
rp = 55500;         % da variare
Vinf = 0.0997;      % da capire
e = 1 + rp*Vinf^2/data.mi;
p = rp*(e + 1);
data.R_SOI = 570000;                                        % [km] from Curtis
data.th0 = -acos((p/data.R_SOI - 1)/e);
data.e_hyp = e;

data.sma = p/(data.e_hyp^2 - 1);                            % [km] sma
data.i = 45*pi/180;                                         % [rad] inclination

Kep0 = [data.sma, data.e_hyp, data.i, 0, 0, data.th0];     % Sc initial condition
[R0, V0] = kep2car(Kep0, data.mi);
M0 = data.M;

data.Y0 = [R0; V0; M0];

data.opt1 = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'Event', @ETfinalEvent);


clearvars -except data Kep0








