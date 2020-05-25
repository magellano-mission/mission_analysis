%% Simulation Setup
data.Norb = 3;                  % number of orbits in Walker constellation
data.walk_phas = 2;             % phase shift in Walker constellation


data.Nsat = 18;                 % number of satellites in Walker constellation
data.bw = 40;                   % beamwidth [deg]
data.inc = deg2rad(45);         % inclination [deg (inserted) --> rad (converted)]
data.sma = 10500;               % semi-major axis [km]


data.Nw = length(data.Nsat);
data.Nb = length(data.bw);
data.Ni = length(data.inc);
data.Na = length(data.sma);

data.trashold = 4; 
data.perturb = false;           % false/true --> perturbations off/on

% Parameters
lon = [-180, 180];              % [deg]  
lat = [-90, 90];                % [deg]
disc = [180, 180];                % discretization grid
data.lon = linspace(lon(1), lon(2), disc(1));
data.lat = linspace(lat(1), lat(2), disc(2));
data.alt = 0;                   % altitude at which evaluate the coverage ( ground level = 0)              
data.N_orbits = 2.7177;              % number of orbits in the simulation

%% Ellipsoid model
data.mi = astroConstants(14);      % Mars planetary constant [km^3/s^2]
data.NT = 1200;                    % timestep

data.rM_eq = 3393.4;               % equatorial radius [km]
rM_pol = 3375.7;                   % polar radius [km]
data.ang_ecc = 2*atan(sqrt((data.rM_eq-rM_pol)/(data.rM_eq+rM_pol)));  % Mars angular eccentricity

T_lla2 = 1.02749125 * 24 * 3600;   % period [s]  
data.wM = 2 * pi / T_lla2;         % [rad/s]
data.theta_Airy_0 = 0;             % Mars principal meridian [rad]

%% plot selection
data.SimType = "GDOP";

clearvars -except data
