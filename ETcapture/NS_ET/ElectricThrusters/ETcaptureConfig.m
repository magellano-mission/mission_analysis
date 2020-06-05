%% Setup of the Electric Thrust capture
% Time
InitData = [2028 06 04 0 0 0 ];           % initial date
data.InitDay = date2mjd2000(InitData);    % initial day

% Propulsion parameter
data.Isp = 4300;                          % specific impulse [s]   
data.g0 = 9.807;                          % [m/s^2]      
data.M = 1600;                            % s/c mass at arrival [kg]

% Planetary parameter
data.mi = astroConstants(14);             % planetary constant [km^3/s^2]
data.switchers.grav = true;               % set gravity effects on/off
data.J2 = 0.00196;                        % Mars J2
data.J3 = 3.1450e-5;                      % Mars J3
data.J4 = -1.5377e-5;                     % Mars J4
T_lla2 = 1.02749125 * 24 * 3600;          % period [s]  
data.wM = 2 * pi / T_lla2;                % Mars angular velocity [rad/s]
V = 1.6318e11;                            % Mars volume [km^3]
data.R_pl = nthroot(3*V/(4*pi), 3);       % Mars equivalent radius [km]
data.R_SOI = 570000;                      % SOI radius [km] (from Curtis)

% Orbit parameters: SOI entry
rp = 59300;                               % pericenter of arrival hyp [km]
Vinf = 0.001;                             % velocity at infinite [km/s]
e = 1 + rp*Vinf^2/data.mi;                % eccentricity of the hyp
p = rp*(e + 1);                           % semilatum rectus [km]
data.th0 = -acos((p/data.R_SOI - 1)/e);   % theta infinite

% Initial conditions
data.e_hyp = e;                           % eccentricity
data.sma = p/(data.e_hyp^2 - 1);          % semi-major axis [km] 
data.i = 45*pi/180;                       % inclination [rad] (from GA)
Kep0 = [data.sma, data.e_hyp, data.i, 0, 0, data.th0];     
[R0, V0] = kep2car(Kep0, data.mi);        % cartesian IC
M0 = data.M;
data.Y0 = [R0; V0; M0];

% Stop conditions
data.thetaStop = 2.2991e+000;                   % [rad] (from GA)
data.timeStop = 4.4092e+006;                 % [s] (from GA)
data.r_final = 12300;                     % NS final radius

% ODE options
data.opt1 = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'Event', @ETParab2EllipseEvent);
data.opt2 = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'Event', @ETCorrectionEvent);
data.opt3 = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
data.opt4 = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'Event', @ETfinalEvent);

clearvars -except data Kep0






