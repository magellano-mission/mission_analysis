%{
script for analyze constellation set-up

%}

%% Matlab Initialization
clear; close all; clc

%% Figure Initialization
set(0,'DefaultFigureUnits', 'normalized');
set(0,'DefaultFigurePosition',[0 0 1 1]);
set(0,'DefaultTextFontSize',18);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultAxesXGrid','on')
set(0,'DefaultAxesYGrid','on')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

%% GNSS orbital parameters

a = 11500;
e = 0;
i = 55;
mi = 42828.3;                   % mars gravity constant [km^2/s^3]
J2 = 1.955e-3;                  % mars J2 gravity coefficient
V = 1.6318e11;                  % mars volume [km^3]
R = nthroot(3*V/(4*pi), 3);     % mars equivalent radius [km]
n = sqrt(mi/a^3);               % mars angular velocity [rad/s]

%% Orbital Parameters secular variation (J2 on GNNS orbit) 
K = -3*J2*sqrt(mi)*R^2/(2*a^(7/2)*(1-e^2)^2);

% variation per seconds in rad
RAAN_dot = K*cosd(i);
PA_dot = K*(5/2*(sind(i))^2 - 2);
theta_dot = n - K*(1 - 3/2*(sind(i))^2);

% variation per days in degrees
RAAN_d = RAAN_dot*86400*180/pi;
PA_d = PA_dot*86400*180/pi;
theta_d = theta_dot*86400*180/pi;

t = 0:365; 

figure; plot(t, t*RAAN_d, 'LineWidth', 2);
xlabel('t [days]'); ylabel(' RAAN [deg]'); 
title('RAAN variation');

figure; plot(t, t*PA_d, 'LineWidth', 2);
xlabel('t [days]'); ylabel(' PA [deg]'); 
title('PA variation');