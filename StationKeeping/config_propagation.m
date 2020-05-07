
%% Simulation Setup
InitData = [2025 01 01 0 0 0 ];                  % initial date
data.InitDay = date2mjd2000(InitData);                   % initial day
data.walker = [1, 1, 1];
% data.bw = 40;              % beamwidth [deg]
data.inc = deg2rad(25);
data.sma = 10500;

% data.trashold = 4; 

% Parameters
% lon = [-180, 180];
% LRG = 15;
% lat = [-LRG, LRG];

% disc = [20 20];                           % discretization grid
% data.lon = linspace(lon(1), lon(2), disc(1));
% data.lat = linspace(lat(1), lat(2), disc(2));
% data.alt = 0;                             % altitude at which evaluate the coverage ( ground level = 0)  

N_orbits = 100;                         % number of orbits
data.mi = astroConstants(14);               % Mars planetary constant
T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
data.tend = 86400*10;

% data.NT = 100000; 
% data.tspan = linspace(0, N_orbits*T_orb, data.NT);

%% perturbation switchers
data.switchers.Moon = false;
data.switchers.grav = true;
data.switchers.SRP = false;

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
Acs = 10;                                           % area exposed to sun
data.m_sat = 250;                                   % satellite mass
data.R_pl = data.rM_eq;                             % planet radius
data.P0 = 4.5*1e-6;                                 % solar radiation pressure coefficient
data.AU = 1.496*1e8;                                % astronomic unit
data.AOverM = Acs/data.m_sat;                       % Am parameter

%% Ephemerides
data.Eph_T0 = date2mjd2000([2024 01 01 0 0 0]);
data.Ephs_Phobos = Ephs;
data.Phobos_time = Time;

Mphobos = 1.07e16; 
Mdeim = 1.4762e15;
G = astroConstants(1);
data.miPhob = Mphobos * G;
data.miDeim = Mdeim * G;

%% ode set
data.opt = odeset('AbsTol', 1e-15, 'RelTol', 1e-13);

clearvars -except data
