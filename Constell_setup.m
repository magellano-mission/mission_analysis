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

a_GNSS = 10500;
e_GNSS = 0;
i_GNSS = 25;
mi = 42828.3;                   % mars gravity constant [km^2/s^3]
J2 = 1.955e-3;                  % mars J2 gravity coefficient
V = 1.6318e11;                  % mars volume [km^3]
R = nthroot(3*V/(4*pi), 3);     % mars equivalent radius [km]
n = sqrt(mi/a_GNSS^3);               % GNSS angular velocity [rad/s]
v_GNSS = sqrt(mi/a_GNSS);            % linear velocity of the GNSS constellation [km/s]

%% Orbital Parameters secular variation (J2 on GNNS orbit) 

a_J2 = 8000:20:11000;

K = -3*J2*sqrt(mi)*R^2./(2*a_J2.^(7/2)*(1-e_GNSS^2)^2);

% variation per seconds in rad
RAAN_dot = K*cosd(i_GNSS);
% PA_dot = K*(5/2*(sind(i_GNSS))^2 - 2);

% variation per days in degrees
RAAN_d = RAAN_dot*86400*365*180/pi;
% PA_d = PA_dot*86400*365*180/pi;

% t = 0:365; 

figure; plot(a_J2, RAAN_d, 'LineWidth', 2);
xlabel('a [km]'); ylabel(' \Delta \Omega [deg/year]'); 
title('RAAN variation due to J2');

% figure; plot(t, t*PA_d, 'LineWidth', 2);
% xlabel('t [days]'); ylabel(' PA [deg]'); 
% title('PA variation');

%% Phasing Maneuvers
% computing the cost of the maneuver wrt the phasing angle and number of revolutions

%%%%%%% NS
T = 2*pi/n;                         % GNSS period [s]
% N = 180;

% dv_phas = NaN(N, N/2);
% for i = 1:N                         % i: phasing angle
%     for j = 3:N/2                     % j: number of revolutions
%         dT_tot = (i*pi/180)/n;
%         dT_rev = dT_tot/j;          
%         T_phas = T - dT_rev;
%         a_phas  = (T_phas*sqrt(mi)/(2*pi))^(2/3);
%         ra_phas = a_GNSS;
%         rp_phas = 2*a_phas - ra_phas;
%         e_phas = (ra_phas - rp_phas)/(ra_phas + rp_phas);
%         v_phas = sqrt(mi/a_phas)*(1 - e_phas);
%         dv_phas(i, j) = 2*abs(v_GNSS - v_phas);
%     end
% end
% 
% surf(1:N, 1:N/2, dv_phas', 'EdgeColor', 'none')
% colorbar; xlabel('phasing angle [deg]'), ylabel('number of revolutions');
% zlabel('\Delta_v [km/s]'), title('phasing maneuver')

m2 = 86400*60;                  % 2 months
nrev = ceil(m2/T);
phas_angle = 60;
dT_tot = (phas_angle*pi/180)/n;
dT_rev = dT_tot/nrev;
T_phas = T - dT_rev;
a_phas  = (T_phas*sqrt(mi)/(2*pi))^(2/3);
ra_phas = a_GNSS;
rp_phas = 2*a_phas - ra_phas;
e_phas = (ra_phas - rp_phas)/(ra_phas + rp_phas);
v_phas = sqrt(mi/a_phas)*(1 - e_phas);
dv_phas_NS = 2*abs(v_GNSS - v_phas)*1e3;


%%%%%%%%% RS
a_RS = 6400;
v_RS = sqrt(mi/a_RS);
phas_angle = 72;
n_RS = sqrt(mi/a_RS^3);
T_RS = 2*pi/n_RS;
nrev = ceil(m2/T_RS);
dT_tot_RS = (phas_angle*pi/180)/n_RS;
dT_rev = dT_tot_RS/nrev;
T_phas = T_RS - dT_rev;
a_phas  = (T_phas*sqrt(mi)/(2*pi))^(2/3);
ra_phas = a_RS;
rp_phas = 2*a_phas - ra_phas;
e_phas = (ra_phas - rp_phas)/(ra_phas + rp_phas);
v_phas = sqrt(mi/a_phas)*(1 - e_phas);
dv_phas_RS = 2*abs(v_RS - v_phas)*1e3;


%% Homhann transfer from RS to ECS
% ECS constellation is not defined, circular orbit in the same plane of RS as an hyphotesys

a_ECS = 6400:50:7400;
N = length(a_ECS);

dv_h = zeros(N, 1);
tof_ECS = zeros(N, 1);
for i = 1:N
    ra_h = a_ECS(i);
    rp_h = a_RS;
    a_h = (ra_h + rp_h)/2;
    v_ECS = sqrt(mi/ra_h);
    vp_h = sqrt(2*mi*(1/rp_h - 1/(2*a_h)));
    va_h = sqrt(2*mi*(1/ra_h - 1/(2*a_h)));
    dv_h(i) = abs(v_ECS - va_h) + abs(v_RS - vp_h);
    tof_ECS(i) = pi*sqrt(a_h^3/mi)/3600;
end

figure; plot(a_ECS, dv_h, 'LineWidth', 2)
xlabel('radius of ECS [km]'), ylabel('\Delta_v [km/s]');
title('homhann to reach ECS')

figure; plot(a_ECS, tof_ECS, 'LineWidth', 2)
xlabel('radius of ECS [km]'), ylabel('\Delta_t [h]');
title('homhann to reach ECS')


%% Homhann transfer from lower orbit to shift RAAN to NS
% ECS constellation is not defined, circular orbit in the same plane of RS as an hyphotesys

a_shift = [8100, 9850];


dv_h = zeros(2, 1);
tof_shift = zeros(2, 1);
for i = 1:2
    ra_h = a_GNSS;
    rp_h = a_shift(i);
    a_h = (ra_h + rp_h)/2;
    v_shift = sqrt(mi/rp_h);
    vp_h = sqrt(2*mi*(1/rp_h - 1/(2*a_h)));
    va_h = sqrt(2*mi*(1/ra_h - 1/(2*a_h)));
    dv_h(i) = abs(v_GNSS - va_h) + abs(v_shift - vp_h);
    tof_shift(i) = pi*sqrt(a_h^3/mi)/3600;
end

figure; plot(a_shift, dv_h, 'o', 'LineWidth', 2)
xlabel('radius of lower orbits [km]'), ylabel('\Delta_v [km/s]');
title('homhann to reach NS')

figure; plot(a_shift, tof_shift, 'o', 'LineWidth', 2)
xlabel('radius of lower orbits [km]'), ylabel('\Delta_t [h]');
title('homhann to reach NS')

%% Low Thrust analytical-approx for ECS transfer from RS orbit
% acceleration along theta direction and constant mass

m_ECS = 1000;                % mass [kg]
T = 0.1:0.005:1;             % thrust [N]
at = T/m_ECS*1e-3;           % acceleration [km/s^2]


figure; hold on

t = 2*(a_ECS(end) - a_RS)./(at*a_RS*sqrt(a_ECS(end)/mi))/86400/30;
plot(T, t, 'LineWidth', 2)

xlabel('on-board continuos thrust [N]'), ylabel('tof [months]');
title('low-thrust to reach ECS')


%% Low Thrust analytical-approx for NS transfer from lower orbit
% acceleration along theta direction and constant mass

m_ECS = 1500;                % mass [kg]
T = 0.01:0.005:1;             % thrust [N]
at = T/m_ECS*1e-3;           % acceleration [km/s^2]

figure; hold on

t1 = 2*(a_GNSS - a_shift(1))./(at*a_shift(1)*sqrt(a_GNSS(end)/mi))/86400/30;
t2 = 2*(a_GNSS - a_shift(2))./(at*a_shift(2)*sqrt(a_GNSS(end)/mi))/86400/30;
plot(T, t1, T, t2, 'LineWidth', 2)

xlabel('on-board continuos thrust [N]'), ylabel('tof [months]');
title('low-thrust to reach ECS')
legend('from a = 8100', 'from a = 9850')


