%% Simulation Setup
% Walker constellation
data.Norb = 3;                  % number of orbits in Walker constellation
data.walk_phas = 2;             % phase shift in Walker constellation

data.Nsat = 21;                 % number of satellites in Walker constellation
data.bw = 40;                   % beamwidth [deg]
data.inc = deg2rad(45);         % inclination [deg (inserted) --> rad (converted)]
data.sma = 12300;               % semi-major axis [km]

data.Nw = length(data.Nsat);
data.Nb = length(data.bw);
data.Ni = length(data.inc);
data.Na = length(data.sma);

% Polar extension
data.PE = false;                        % add or not the polar extention
if data.PE == true
    data.Nsat_pol = 6;
    data.bw_pol = data.bw;              % beamwidth [deg]
    data.inc_pol = deg2rad(90);
    data.sma_pol = data.sma + 100;      % still have to understand how to avoid 
end                                     % sats collisions between wlk and polar

data.mask = 12;                  % minimum elevation angle [deg] for ground users
% data.mask = 0;                     % minimum elevation angle [deg] for satellites
data.trashold = 4; 
data.perturb = false;              % false/true --> perturbations off/on

% Parameters
lon = [-180, 180];                 % [deg]  
lat = [-90, 90];                   % [deg]

data.alt = 0;                    % altitude at which evaluate the coverage ( ground level = 0)              
% data.alt = 1.5150e+003;          % RS altitude
% data.alt = 4.0150e+003;
data.N_orbits = 2.1435;            % number of orbits in the simulation

% data.study = "DOP";
data.study = "coverage";

if data.study == "DOP"
    disc = [160, 150];             % discretization grid
else
    disc = [90, 90]; 
end

data.lon = linspace(lon(1), lon(2), disc(1));
data.lat = linspace(lat(1), lat(2), disc(2));

%% Ellipsoid model
data.mi = astroConstants(14);      % Mars planetary constant [km^3/s^2]
data.NT = 900;                    % timestep

data.rM_eq = 3393.4;               % equatorial radius [km]
rM_pol = 3375.7;                   % polar radius [km]
data.ang_ecc = 2*atan(sqrt((data.rM_eq-rM_pol)/(data.rM_eq+rM_pol)));  % Mars angular eccentricity

T_lla2 = 1.02749125 * 24 * 3600;   % period [s]  
data.wM = 2 * pi / T_lla2;         % [rad/s]
data.theta_Airy_0 = 0;             % Mars principal meridian [rad]

%% plot selection
if data.study == "DOP"
    data.SimType = "GDOP";
else
    % data.SimType = "Dani";
    data.SimType = "FPCP"; 
    %      data.SubType = "mean"; 
         data.SubType = "min"; 
    % data.SimType = 'none';
    % data.SimType = "percent_cov";
    % data.SimType = "max_time_gap";
    % data.SimType = "mean_time_gap";
    % data.SimType = "plot_belli";
    % data.SimType =  "Time_varying";
    
end




clearvars -except data
