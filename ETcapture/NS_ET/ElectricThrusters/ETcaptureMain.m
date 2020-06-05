% ET Capture
close all; clear; clc

%% Adding to path 
addpath(genpath(fileparts(strcat(pwd, '/functions'))))

%% Figure Initialization
set(0,'DefaultFigureUnits', 'normalized');
set(0,'DefaultFigurePosition', [0 0 1 1]);
set(0,'DefaultTextFontSize', 18);
set(0,'DefaultAxesFontSize', 18);
set(0,'DefaultAxesXGrid', 'on')
set(0,'DefaultAxesYGrid', 'on')
set(0,'defaultLegendInterpreter', 'latex');
set(0,'defaultAxesTickLabelInterpreter', 'latex');

%% Importing Data
run('ETcaptureConfig.m')                   

%% Computations
tic
TT = 56.8800e+006;                     % long period, Event function inside the ode
data.T = -0.196;                % [N] (for clearness the Thrust is here)
data.flag = 0;

% NS capture maneuver
data.direction = "tangential";
[T, Y] = ode113(@ETcaptIntegration, [0, TT], data.Y0, data.opt4, data);
parameters = car2kep(Y(end,1:3),Y(end,4:6),data.mi)
r_fin = norm(Y(end,1:3))
mass = Y(end,7)

%% Post-process
vecDates = data.InitDay +  T/86400;

vecSun = zeros(length(vecDates),3);
for kk = 1:length(vecDates)
   vecSun1 = uplanet(vecDates(kk),4); 
   [vecSun(kk,:), ~] = kep2car(vecSun1,data.mi);
end

incMars = 1.85061*pi/180;   % [rad]
A = [cos(incMars) sin(incMars) 0
     -sin(incMars) cos(incMars) 0
     0 0 1];

Sun2Mars = zeros(length(T),3);
for jj = 1:length(T)
    Sun2Mars(jj,:) = A*vecSun(jj,1:3)';
end

Mars2sc = Y(:,1:3);

SunDirection = - (Mars2sc + Sun2Mars);


%% Graphs
% Orbit
figure()
rM = almanac('mars','radius','kilometers','sphere');
I = imread('Mars.jpg');                            
RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];                       % Mars image x sizes
RI.YWorldLimits = [-90 90];                         % Mars image y sizes
[XM, YM, ZM] = ellipsoid(0, 0, 0, rM, rM, rM, 100); % spheric centered Mars
planet = surf(XM, YM, -ZM,'Edgecolor', 'none');
hold on
set(planet,'FaceColor','texturemap','Cdata',I)
set(get(get(planet,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(gca,'Color','black')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')

plot3(Y(:, 1), Y(:, 2), Y(:, 3),'Color',[155/255, 155/255, 155/255]), hold on, axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
% [0.9490, 0.4745, 0.3137]
% 
% [0.1020, 0.6667, 0.74120]
% 
% [155/255, 155/255, 155/255]

% other stuff for check

Yk = zeros(length(Y),6);
for kk = 1:length(T)
    Yk(kk,1:6) = car2kep(Y(kk,1:3),Y(kk,4:6),data.mi);
end

figure()
subplot(2,3,1)
plot(T,Yk(:,1))
ylabel('SMA')

subplot(2,3,2)
plot(T,Yk(:,2))
ylabel('ecc')

subplot(2,3,3)
plot(T,Yk(:,3)*180/pi)
ylabel('incl')

subplot(2,3,4)
plot(T,Yk(:,4)*180/pi)
ylabel('RAAN')

subplot(2,3,5)
plot(T,Yk(:,5)*180/pi)
ylabel('omeghino')

subplot(2,3,6)
plot(T,Yk(:,6)*180/pi)
ylabel('theta')



