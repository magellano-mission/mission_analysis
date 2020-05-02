
%% Ellipsoid model
para.rM_eq = 3393.4;         % equatorial radius [km]
rM_pol = 3375.7;                   % polar radius [km]
para.ang_ecc = 2*atan(sqrt((para.rM_eq-rM_pol)/(para.rM_eq+rM_pol)));  % Mars angular eccentricity

T_lla2 = 1.02749125 * 24 * 3600;   % period [s]  
para.wM = 2 * pi / T_lla2;              % [rad/s]
para.theta_Airy_0 = 0;             % Mars principal meridian

%% Setup
wlk_vec = [15];       % number of satellites
bw = 20;           % beamwidth [deg]
inclinations = deg2rad(55);
semi_major_axes = 11500;
Treshold = 4; 

% Parameters
lon = [-180, 180];          
lat = [-90, 90];
disc = [15 6];             % discretization grid
lon = linspace(lon(1), lon(2), disc(1));
lat = linspace(lat(1), lat(2), disc(2));
alt = 0;                   % altitude at which evaluate the coverage ( ground level = 0)
timesteps = 1000;               
N_orbits = 3;              % number of orbits

%% plot selection

% SimType = "Dani";
SimType = "FPCP";
