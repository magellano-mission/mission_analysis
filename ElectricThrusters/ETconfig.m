%%
InitData = [2025 01 01 0 0 0 ];                         % initial date
data.InitDay = date2mjd2000(InitData);                  % initial day

data.mi = astroConstants(14);
V = 1.6318e11;                  % mars volume [km^3]
data.R_pl = nthroot(3*V/(4*pi), 3);     % mars equivalent radius [km]

data.M = 250;            % kg

data.dth = pi;
data.sma = 6400;
data.i = 0*pi/180;
Kep0 = [data.sma, 0, data.i , 0, 0, 0];
[R0, V0] = kep2car(Kep0, data.mi);
data.Y0 = [R0; V0];

Kepfree = [data.sma, 0, data.i, 0, 0, data.dth];
[R0, V0] = kep2car(Kepfree, data.mi);
data.Yfree = [R0; V0];

data.opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

clearvars -except data Kep0