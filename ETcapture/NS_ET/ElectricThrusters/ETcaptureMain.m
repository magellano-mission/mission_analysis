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
TT = 56.9339e+006;                     % long period, Event function inside the ode
data.T = -0.2;                % [N] (for clearness the Thrust is here)
data.flag = 0;

%% Try for the GA and set of the initial conditions
% data.direction = "tangential";
% [T, Y] = ode113(@ETcaptIntegration, [0, TT], data.Y0, data.opt4, data);
% parameters = car2kep(Y(end,1:3),Y(end,4:6),data.mi)
% r_fin = norm(Y(end,1:3))
% mass = Y(end,7)

%% Integration until ellipse
data.direction = "tangential";
[T1, Y1] = ode113(@ETcaptIntegration, [0, TT], data.Y0, data.opt1, data);

%% Integration until reach descending or ascending node
data.T = 0;           % coasting maneuver, thrust is off
kepDesc = car2kep(Y1(end,1:3),Y1(end,4:6), data.mi);
kepDesc(7) = Y1(end,7);
[T2, Y2] = ode113(@GaussIntegration, [T1(end), TT], kepDesc, data.opt2, data);

%% Integration for the change of RAAN
data.direction = "perpendicular";
data.T = -0.2; 
kepChangeRAAN = Y2(end,:);
kepChangeRAAN(7) = Y2(end,7);
[T3, Y3] = ode113(@GaussIntegration, [T2(end), T2(end) + data.timeStop], kepChangeRAAN, data.opt3, data);

%% Integration until close the orbit
data.direction = "tangential";
data.T = -0.2;
[RR, VV] = kep2car(Y3(end,1:6),data.mi);
cartIC = [RR; VV; Y3(end,7)];
[T4, Y4] = ode113(@ETcaptIntegration, [T3(end), TT], cartIC, data.opt4, data);

%% Post-process
% r_fin = norm(Y1(end,1:3))

Y2kep = Y2(:,1:6);
Y2cart = zeros(length(T2),3);
for kk = 1: length(T2)
    Y2cart(kk,:) = kep2car(Y2kep(kk,:), data.mi);
end

Y3kep = Y3(:,1:6);
Y3cart = zeros(length(T3),3);
for kk = 1: length(T3)
    Y3cart(kk,:) = kep2car(Y3kep(kk,:), data.mi);
end


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
%set(gca,'Color','black')
plot3(Y1(:, 1), Y1(:, 2), Y1(:, 3),'Colo',[155/255, 155/255, 155/255]), hold on, axis equal
plot3(Y2cart(:, 1), Y2cart(:, 2), Y2cart(:, 3),'k')
plot3(Y3cart(:, 1), Y3cart(:, 2), Y3cart(:, 3),'Color',[0.1020, 0.6667, 0.74120])
plot3(Y4(:, 1), Y4(:, 2), Y4(:, 3),'Color', [0.9490,0.4745,0.3137])
% set(gca,'Visible','off')

% [0.9490, 0.4745, 0.3137]
% 
% [0.1020, 0.6667, 0.74120]
% 
% [155/255, 155/255, 155/255]
% Mass decrease
% hold off
% figure()
% plot(T_tot,mass)
    
%% EXTRA for DAG - da sistemare nel caso di utilizzo

% data.flag = 1;        this is only in case of on/off DAG

% % Integration for the change of RAAN
% data.flag = 1;
% data.Y02 = Y6(end,:);
% data.kep02 = car2kep(data.Y02(1:3),data.Y02(4:6), data.mi);
% [T7, Y7] = ode113(@ETcaptIntegration, [T6(end), TT], data.Y02, data.opt5, data);

Y4k = zeros(length(Y4),6);
for kk = 1:length(T4)
    Y4k(kk,1:6) = car2kep(Y4(kk,1:3),Y4(kk,4:6),data.mi);
end

figure
subplot(2,3,1)
plot(T4,Y4k(:,1))
ylabel('SMA')

subplot(2,3,2)
plot(T4,Y4k(:,2))
ylabel('ecc')

subplot(2,3,3)
plot(T4,Y4k(:,3)*180/pi)
ylabel('incl')

subplot(2,3,4)
plot(T4,Y4k(:,4)*180/pi)
ylabel('RAAN')

subplot(2,3,5)
plot(T4,Y4k(:,5)*180/pi)
ylabel('omeghino')

subplot(2,3,6)
plot(T4,Y4k(:,6)*180/pi)
ylabel('theta')



