%% Simulation Setup
data.Norb = 1;                  %wl
data.walk_phas = 1;


data.Nsat = 6;
data.bw = 97;              % beamwidth [deg]
data.inc = deg2rad(0);
data.sma = 4430;


data.Nw = length(data.Nsat);
data.Nb = length(data.bw);
data.Ni = length(data.inc);
data.Na = length(data.sma);

data.trashold = 1; 
data.perturb = false;

% Parameters
lon = [-180, 180];          
lat = [-40, 40];
disc = [70, 80];             % discretization grid
data.lon = linspace(lon(1), lon(2), disc(1));
data.lat = linspace(lat(1), lat(2), disc(2));
data.alt = 0;                  % altitude at which evaluate the coverage ( ground level = 0)              
data.N_orbits = 2;             % number of orbits

%% Ellipsoid model
data.mi = astroConstants(14);            % Mars planetary constant
data.NT = 700; 

data.rM_eq = 3393.4;         % equatorial radius [km]
rM_pol = 3375.7;             % polar radius [km]
data.ang_ecc = 2*atan(sqrt((data.rM_eq-rM_pol)/(data.rM_eq+rM_pol)));  % Mars angular eccentricity

T_lla2 = 1.02749125 * 24 * 3600;   % period [s]  
data.wM = 2 * pi / T_lla2;         % [rad/s]
data.theta_Airy_0 = 0;             % Mars principal meridian

%% plot selection

% data.SimType = "Dani";
% data.SimType = "FPCP"; 
%      data.SubType = "mean"; 
%      data.SubType = "min"; 
% data.SimType = 'none';
% data.SimType = "percent_cov";
% data.SimType = "max_time_gap";
% data.SimType = "mean_time_gap";
% data.SimType = "plot_belli";
 data.SimType =  "Time_varying";

clearvars -except data
