%% Setup of the Electric Thrust capture
% Time
InitData = [2028 06 04 0 0 0 ];           % initial date
data.InitDay = date2mjd2000(InitData);    % initial day

% Propulsion parameter
data.Isp_fun = @(X) 0.051346711071873*X.^2-6.382637268333952*X+3554.166890279770;
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
data.Vinf = 0.001;                        % velocity at infinite [km/s]
data.i = 45*pi/180;                       % inclination [rad] 
data.Om0 = 359.999*pi/180;

% Stop conditions
data.r_final = 12300;                     % NS final radius
data.eccF = 0.0005;

% ODE options
data.opt = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'Event', @ETfinalEvent);

clearvars -except data Kep0






