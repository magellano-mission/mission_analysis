%% Setup of the Electric Thrust capture
% Time
SM = "NS1";
% SM = "NS2";
% SM = "NS3";
% SM = "ECS";
% SM = "RSeq";
% SM = "RSpol";

% Propulsion parameter
data.Isp_fun = @(X) 0.051346711071873*X.^2-6.382637268333952*X+3554.166890279770;
data.g0 = 9.807;                          % [m/s^2]  

switch SM
    case "NS1"
        data.M = 1750;                            % s/c mass at arrival [kg]
        data.Apanels = 52;
        data.r_final = 12300;                     % NS final radius
        data.i = 45*pi/180;                       % inclination [rad]
        InitData = [2028 06 12 0 0 0 ];           % initial date
        data.InitDay = date2mjd2000(InitData);    % initial day
        data.rp0 = [50000, 62000];
        data.rp = 58628;
    case "ECS"
        data.M = 1450;
        data.Apanels = 52;
        data.r_final = 7400;                      % ECS final radius
        data.i = 0*pi/180;                        % inclination [rad]
        InitData = [2028 06 16 0 0 0 ];           % initial date
        data.InitDay = date2mjd2000(InitData);    % initial day
        data.rp0 = [32000, 45000];
        data.rp = 45562;
    case "RSeq"
        data.M = 1350;
        data.Apanels = 52;
        data.r_final = 4900;                      % NS final radius
        data.i = 0*pi/180;                        % inclination [rad]
        InitData = [2028 06 14 0 0 0 ];           % initial date
        data.InitDay = date2mjd2000(InitData);    % initial day
        data.rp0 = [28000, 40000];
        data.rp = 37654.4;
    case "NS2"
        data.M = 1750;                            % s/c mass at arrival [kg]
        data.Apanels = 52;
        data.r_final = 12300;                     % NS final radius
        data.i = 45*pi/180;                       % inclination [rad]
        InitData = [2028 12 13 0 0 0 ];           % initial date
        data.InitDay = date2mjd2000(InitData);    % initial day
        data.rp0 = [50000, 65000];
        data.rp = 54657;
    case "RSpol"
        data.M = 1350;
        data.Apanels = 52;
        data.r_final = 4900;                      % NS final radius
        data.i = 89.999*pi/180;                        % inclination [rad]
        InitData = [2028 12 16 0 0 0 ];           % initial date
        data.InitDay = date2mjd2000(InitData);    % initial day
        data.rp0 = [30000, 45000];
        data.rp = 35119;
    case "NS3"
        data.M = 1750;                            % s/c mass at arrival [kg]
        data.Apanels = 52;
        data.r_final = 12300;                     % NS final radius
        data.i = 45*pi/180;                       % inclination [rad]
        InitData = [2029 06 22 0 0 0 ];           % initial date
        data.InitDay = date2mjd2000(InitData);    % initial day 
        data.rp0 = [50000, 62000];
        data.rp = 56723;
end
        
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
data.Om0 = 359.999*pi/180;

% ODE options
data.opt = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'Event', @ETfinalEvent);

clearvars -except data Kep0






