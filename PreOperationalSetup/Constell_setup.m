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

a_GNSS = 12300;
e_GNSS = 0;
i_GNSS = 45;
mi = 42828.3;                       % mars gravity constant [km^2/s^3]
J2 = 1.955e-3;                      % mars J2 gravity coefficient
V = 1.6318e11;                      % mars volume [km^3]
R = nthroot(3*V/(4*pi), 3);         % mars equivalent radius [km]
n = sqrt(mi/a_GNSS^3);              % GNSS angular velocity [rad/s]
v_GNSS = sqrt(mi/a_GNSS);           % linear velocity of the GNSS constellation [km/s]


%% Phasing Maneuvers
% computing the cost of the maneuver wrt the phasing angle and number of revolutions

%%%%%%% NS
T = 2*pi/n;                         % GNSS period [s]

m2 = 86400*90;                      % 3 months
nrev = ceil(m2/T);
phas_angle = 51.43;
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
a_RS = 4900;
v_RS = sqrt(mi/a_RS);
phas_angle = 90;
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


%%%%%%%%% ECS
a_ECS = 7400;
v_ECS = sqrt(mi/a_ECS);
phas_angle = 180;
n_ECS = sqrt(mi/a_ECS^3);
T_ECS = 2*pi/n_ECS;
nrev = ceil(m2/T_ECS);
dT_tot_ECS = (phas_angle*pi/180)/n_ECS;
dT_rev = dT_tot_ECS/nrev;
T_phas = T_ECS - dT_rev;
a_phas  = (T_phas*sqrt(mi)/(2*pi))^(2/3);
ra_phas = a_ECS;
rp_phas = 2*a_phas - ra_phas;
e_phas = (ra_phas - rp_phas)/(ra_phas + rp_phas);
v_phas = sqrt(mi/a_phas)*(1 - e_phas);
dv_phas_ECS = 2*abs(v_ECS - v_phas)*1e3;


%% phasing plot PM3

I = imread('mars.jpg');                             % earth image
RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];                       % earth image x sizes
RI.YWorldLimits = [-90 90];                         % earth image y sizes

figure; hold on; view(0, 90)
set(gca,'Color','none','visible','off','YDir','normal')
axis equal
% axis([-1.3e4, 1.3e4, -1.3e4, 1.3e4, -1.3e4, 1.3e4])

[X, Y, Z] = ellipsoid(0, 0, 0, R, R, R, 100); % spheric centered earth
planet = surf(X, Y, -Z,'Edgecolor', 'none');
set(planet,'FaceColor','texturemap','Cdata',I)

OrbPar1 = [a_GNSS, 0, 25*180/pi, 0, 0, 0];

[~, R1] = PlotConic(OrbPar1, [150 150 150]/255, 1, '-', 1:360);
plot3(R1(1, 1), R1(1, 2), R1(1, 3), 'o', 'MarkerEdgeColor', [0.9490    0.4745    0.3137], 'MarkerFaceColor',...
    [0.9490    0.4745    0.3137], 'MarkerSize', 20 )
plot3(R1(61:60:end, 1), R1(61:60:end, 2), R1(61:60:end, 3), 'o', 'MarkerEdgeColor', [0.1020    0.6667    0.74120], 'MarkerFaceColor',...
    [0.1020    0.6667    0.74120], 'MarkerSize', 7)

OrbPar1 = [a_GNSS, 0, 25*180/pi, 120*pi/180, 0, 0];

[~, R1] = PlotConic(OrbPar1, [150 150 150]/255, 1, '-', 1:360);
plot3(R1(1, 1), R1(1, 2), R1(1, 3), 'o', 'MarkerEdgeColor', [0.9490    0.4745    0.3137], 'MarkerFaceColor',...
    [0.9490    0.4745    0.3137], 'MarkerSize', 20 )
plot3(R1(61:60:end, 1), R1(61:60:end, 2), R1(61:60:end, 3), 'o', 'MarkerEdgeColor', [0.1020    0.6667    0.74120], 'MarkerFaceColor',...
    [0.1020    0.6667    0.74120], 'MarkerSize', 7)

OrbPar1 = [a_GNSS, 0, 25*180/pi, 240*pi/180, 0, 0];

[~, R1] = PlotConic(OrbPar1, [150 150 150]/255, 1, '-', 1:360);
plot3(R1(1, 1), R1(1, 2), R1(1, 3), 'o', 'MarkerEdgeColor', [0.9490    0.4745    0.3137], 'MarkerFaceColor',...
    [0.9490    0.4745    0.3137], 'MarkerSize', 20 )
plot3(R1(61:60:end, 1), R1(61:60:end, 2), R1(61:60:end, 3), 'o', 'MarkerEdgeColor', [0.1020    0.6667    0.74120], 'MarkerFaceColor',...
    [0.1020    0.6667    0.74120], 'MarkerSize', 7)

%% Homhann transfer from lower orbit to shift RAAN to NS
% ECS constellation is not defined, circular orbit in the same plane of RS as an hyphotesys

dv_h = zeros(2, 1);
tof_shift = zeros(2, 1);
for i = 1:2
    ra_h = a_GNSS;
    rp_h = a_shift(i);
    a_h(i) = (ra_h + rp_h)/2;
    e_h(i) = abs((ra_h - rp_h)/(ra_h + rp_h));
    v_shift = sqrt(mi/rp_h);
    vp_h = sqrt(2*mi*(1/rp_h - 1/(2*a_h(i))));
    va_h = sqrt(2*mi*(1/ra_h - 1/(2*a_h(i))));
    dv_h(i) = abs(v_GNSS - va_h) + abs(v_shift - vp_h);
    tof_shift(i) = pi*sqrt(a_h(i)^3/mi)/3600;
end

I = imread('mars.jpg');                             % earth image
RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];                       % earth image x sizes
RI.YWorldLimits = [-90 90];                         % earth image y sizes

figure; hold on; view(-140, 30)
set(gca,'Color','none','visible','off','YDir','normal')
axis equal
% axis([-1.3e4, 1.3e4, -1.3e4, 1.3e4, -1.3e4, 1.3e4])

[X, Y, Z] = ellipsoid(0, 0, 0, R, R, R, 100); % spheric centered earth
planet = surf(X, Y, -Z,'Edgecolor', 'none');
set(planet,'FaceColor','texturemap','Cdata',I)

OrbPar1 = [a_shift(1), 0, 25*pi/180, 240*pi/180, 0, 0];
OrbParh1 = [a_h(1), e_h(1), 25*pi/180, 240*pi/180, 0, 0];
OrbPar2 = [a_shift(2), 0, 25*pi/180, 120*pi/180, 0, 0];
OrbParh2 = [a_h(2), e_h(2), 25*pi/180, 120*pi/180, 0, 0];
OrbParf1 = [a_GNSS, 0, 25*pi/180, 240*pi/180, 0, 0];
OrbParf2 = [a_GNSS, 0, 25*pi/180, 120*pi/180, 0, 0];

[h1, R1] = PlotConic(OrbPar1, [0.9490    0.4745    0.3137], 1, '--', 1:360);
[h2, R2] = PlotConic(OrbPar2, [0.1020    0.6667    0.74120], 1, '--', 1:360);
[h_h1, Rh1] = PlotConic(OrbParh1, [0.9490    0.4745    0.3137], 2, '-', 1:180);
[h_h2, Rh2] = PlotConic(OrbParh2, [0.1020    0.6667    0.74120], 2, '-', 1:180);
[hf1, Rf1] = PlotConic(OrbParf1, [0.9490    0.4745    0.3137], 1, '--', 1:360);
[hf2, Rf2] = PlotConic(OrbParf2, [0.1020    0.6667    0.74120], 1, '--', 1:360);


