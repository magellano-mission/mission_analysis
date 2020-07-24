
%% Simulation Setup
InitData = [2025 01 01 0 0 0 ];                         % initial date
data.InitDay = date2mjd2000(InitData);                  % initial day
data.walker = [1, 1, 1];
data.inc = deg2rad(1e-3);
data.sma = 4900;

data.X0_kep = [data.sma, 1e-5, data.inc, 89.9/180*pi, pi/3, 0]; 

data.mi = astroConstants(14);                           % Mars planetary constant
data.tend = 86400*365*50;

%% perturbation switchers
% data.switchers.Moon = true;
data.switchers.grav = true;
data.switchers.SRP = true;

%% Mars model
data.rM_eq = 3393.4;         % equatorial radius [km]
rM_pol = 3375.7;                   % polar radius [km]
data.ang_ecc = 2*atan(sqrt((data.rM_eq-rM_pol)/(data.rM_eq+rM_pol)));  % Mars angular eccentricity

T_lla2 = 1.02749125 * 24 * 3600;   % period [s]  
data.wM = 2 * pi / T_lla2;         % [rad/s]
data.theta_Airy_0 = 0;             % Mars principal meridian

data.J2 = 0.00196; %Mars J2
data.J3 = 3.1450e-5; %Mars J3
data.J4 = -1.5377e-5; %Mars J4

%% SRP
data.CR = 1.5;                                      % reflection coefficient
Acs = 15;                                           % area exposed to sun
data.m_sat = 180;                                   % satellite mass
data.R_pl = data.rM_eq;                             % planet radius
data.S0 = 63.15e6*696000^2/2.998e8;                 % solar radiation pressure coefficient
data.eps = 5.65;                                    % Mars ecliptic plane[deg]
data.AOverM = Acs/data.m_sat;                       % Am parameter

% %% Ephemerides
% data.Eph_T0 = date2mjd2000([2024 01 01 0 0 0]);
% Eph_time_mjd2000 = Time/86400 + data.Eph_T0;
% mjd2000_end = data.InitDay + data.tend/86400;
% 
% if (mjd2000_end > Eph_time_mjd2000(end)) && data.switchers.Moon
%     error('simulation time is too long, not compatible with Moon Ephemerides')
% end
% 
% ind_start = find(Eph_time_mjd2000 > data.InitDay);
% ind_start = ind_start(1) - 1;
% 
% ind_end = find(Eph_time_mjd2000 > mjd2000_end);
% ind_end = ind_end(1) + 1;
% 
% data.Ephs_Phobos = Ephs(ind_start:ind_end, :);
% data.Phobos_time_mjd2000 = Eph_time_mjd2000(ind_start:ind_end);
% 
% Mphobos = 1.07e16; 
% Mdeim = 1.4762e15;
% G = astroConstants(1);
% data.miPhob = Mphobos * G;
% data.miDeim = Mdeim * G;

%% ode set
data.opt = odeset('AbsTol', 1e-9, 'RelTol', 1e-7);

clearvars -except data
