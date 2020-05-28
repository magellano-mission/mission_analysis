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

%% Integration until ellipse
data.direction = "tangential";
[T1, Y1] = ode113(@ETcaptIntegration, [0, TT], data.Y0, data.opt4, data);

%% Post-process
r_fin = norm(Y1(end,1:3))

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
plot3(Y1(:, 1), Y1(:, 2), Y1(:, 3),'Color',[0.9490,0.4745,0.3137]), axis equal, hold on
set(gca,'Visible','Off')

% Mass decrease
hold off
figure()
plot(T1,Y1(:,7))

Y1k = zeros(length(Y1),6);
for kk = 1:length(T1)
    Y1k(kk,1:6) = car2kep(Y1(kk,1:3),Y1(kk,4:6),data.mi);
end

figure
subplot(2,3,1)
plot(T1,Y1k(:,1))
ylabel('SMA')

subplot(2,3,2)
plot(T1,Y1k(:,2))
ylabel('ecc')

subplot(2,3,3)
plot(T1,Y1k(:,3))
ylabel('incl')

subplot(2,3,4)
plot(T1,Y1k(:,4))
ylabel('RAAN')

subplot(2,3,5)
plot(T1,Y1k(:,5))
ylabel('omeghino')

subplot(2,3,6)
plot(T1,Y1k(:,6))
ylabel('theta')
    
%% EXTRA for DAG - da sistemare nel caso di utilizzo

% data.flag = 1;        this is only in case of on/off DAG

% % Integration for the change of RAAN
% data.flag = 1;
% data.Y02 = Y6(end,:);
% data.kep02 = car2kep(data.Y02(1:3),data.Y02(4:6), data.mi);
% [T7, Y7] = ode113(@ETcaptIntegration, [T6(end), TT], data.Y02, data.opt5, data);





