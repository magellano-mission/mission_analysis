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

a_GNSS = 11500;
e_GNSS = 0;
i_GNSS = 55;
mi = 42828.3;                   % mars gravity constant [km^2/s^3]
J2 = 1.955e-3;                  % mars J2 gravity coefficient
V = 1.6318e11;                  % mars volume [km^3]
R = nthroot(3*V/(4*pi), 3);     % mars equivalent radius [km]
n = sqrt(mi/a_GNSS^3);               % mars angular velocity [rad/s]
v_GNSS = sqrt(mi/a_GNSS);            % linear velocity of the GNSS constellation [km/s]

%% Orbital Parameters secular variation (J2 on GNNS orbit) 
K = -3*J2*sqrt(mi)*R^2/(2*a_GNSS^(7/2)*(1-e_GNSS^2)^2);

% variation per seconds in rad
RAAN_dot = K*cosd(i_GNSS);
PA_dot = K*(5/2*(sind(i_GNSS))^2 - 2);

% variation per days in degrees
RAAN_d = RAAN_dot*86400*180/pi;
PA_d = PA_dot*86400*180/pi;

t = 0:365; 

figure; plot(t, t*RAAN_d, 'LineWidth', 2);
xlabel('t [days]'); ylabel(' RAAN [deg]'); 
title('RAAN variation');

figure; plot(t, t*PA_d, 'LineWidth', 2);
xlabel('t [days]'); ylabel(' PA [deg]'); 
title('PA variation');

%% Phasing Maneuver in GNSS orbit
% computing the cost of the maneuver wrt the phasing angle and number of revolutions

T = 2*pi/n;                         % GNSS period [s]
N = 180;

dv_phas = NaN(N, N/2);
for i = 1:N                         % i: phasing angle
    for j = 3:N/2                     % j: number of revolutions
        dT_tot = (i*pi/180)/n;
        dT_rev = dT_tot/j;          
        T_phas = T - dT_rev;
        a_phas  = (T_phas*sqrt(mi)/(2*pi))^(2/3);
        ra_phas = a_GNSS;
        rp_phas = 2*a_phas - ra_phas;
        e_phas = (ra_phas - rp_phas)/(ra_phas + rp_phas);
        v_phas = sqrt(mi/a_phas)*(1 - e_phas);
        dv_phas(i, j) = 2*abs(v_GNSS - v_phas);
    end
end

surf(1:N, 1:N/2, dv_phas', 'EdgeColor', 'none')
colorbar; xlabel('phasing angle [deg]'), ylabel('number of revolutions');
zlabel('\Delta_v [km/s]'), title('phasing maneuver')

%% Homhann transfer from ECS to GNSS
% ECS constellation is not defined, circular orbit in the same plane of GNSS as an hyphotesys

a_ECS = a_GNSS:500:25000;
N = length(a_ECS);

dv_h = zeros(N, 1);
for i = 1:N
    ra_h = a_ECS(i);
    rp_h = a_GNSS;
    a_h = (ra_h + rp_h)/2;
    v_ECS = sqrt(mi/a_ECS(i));
    vp_h = sqrt(2*mi*(1/rp_h - 1/(2*a_h)));
    va_h = sqrt(2*mi*(1/ra_h - 1/(2*a_h)));
    dv_h(i) = abs(v_ECS - va_h) + abs(v_GNSS - vp_h);
end

figure; plot(a_ECS, dv_h, 'LineWidth', 2)
xlabel('radius of ECS orbit [km]'), ylabel('\Delta_v [km/s]');
title('ECS to GNSS transfer')

%% Low Thrust analytical-approx for GNSS to ECS transfer
% acceleration along theta direction and constant mass

m_GMS = 1000;               % mass [kg]
T = 0.01:0.005:1;           % thrust [N]
a_GMS = R + 200;            % radius of the final orbit [km]
at = T/m_GMS*1e-3;          % acceleration [km/s^2]

t = (a_GNSS - a_GMS)./(at*a_GMS*sqrt(a_GNSS/mi))/86400/30;

figure; plot(T, t, 'LineWidth', 2)
xlabel('on-board continuos thrust [N]'), ylabel('tof [months]');
title('GNSS to GMS transfer')

% Tlt = 2*12*30*86400;                    % transfer time wanted [s]
% 
% a_GMS = a_GNSS/(1 - at*Tlt*sqrt(a_GNSS/mi));


