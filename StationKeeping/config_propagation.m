
%% Simulation Setup
data.walker = [15, 3, 2];
data.bw = 40;              % beamwidth [deg]
data.inc = deg2rad(25);
data.sma = 10500;

data.trashold = 4; 

% Parameters
lon = [-180, 180];
LRG = 15;
lat = [-LRG, LRG];

disc = [20 20];             % discretization grid
data.lon = linspace(lon(1), lon(2), disc(1));
data.lat = linspace(lat(1), lat(2), disc(2));
data.alt = 0;                   % altitude at which evaluate the coverage ( ground level = 0)  

% variations
data.adot = 0;
data.edot = 0;
data.idot = 1e-7;

data.N_orbits = 100;            % number of orbits
data.mi = astroConstants(14);            % Mars planetary constant
data.T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
data.NT = 10000; 
data.tspan = linspace(0, data.N_orbits*data.T_orb, data.NT);
data.n = sqrt(data.mi/data.sma^3);             % angular velocity

%% Ellipsoid model
data.rM_eq = 3393.4;         % equatorial radius [km]
rM_pol = 3375.7;                   % polar radius [km]
data.ang_ecc = 2*atan(sqrt((data.rM_eq-rM_pol)/(data.rM_eq+rM_pol)));  % Mars angular eccentricity

T_lla2 = 1.02749125 * 24 * 3600;   % period [s]  
data.wM = 2 * pi / T_lla2;         % [rad/s]
data.theta_Airy_0 = 0;             % Mars principal meridian

J2 = 0.00196; %Mars J2
J3 = 3.1450e-5; %Mars J3
J4 = -1.5377e-5; %Mars J4

%% Ephemerides
% data.Eph_T0 = mjd2000([2024 - 01 -01])

clearvars -except data
