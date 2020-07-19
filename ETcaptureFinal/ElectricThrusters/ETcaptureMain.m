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
set(0,'defaultLegendInterpreter', 'none');
set(0,'defaultAxesTickLabelInterpreter', 'none');

%% Importing Data
run('ETcaptureConfig.m')                   

%% Computations
tic
data.TT = 1e8;              % long period, Event function inside the ode
data.direction = "tangential";

%% Optimization --> results reported in ETcaptureConfig.m
options = optimset('Display','iter','TolX',1e-4);
rp = fzero(@(x) BCfindRHyp(x, data), data.rp0, options);

%% Retrieve trajectory
[~, Y, T, VectorThrust, VectorLight, VectorThRange] = BCfindRHyp(rp, data);

%% Post-process
parameters = car2kep(Y(end,1:3),Y(end,4:6),data.mi);
mass = Y(end,7);
TOF = T(end)/86400;
ArrivalDate = mjd20002date(data.InitDay + TOF);

timeInterval = T(2:end)-T(1:end-1);
Itot = sum(abs(VectorThrust(1:end-1)).*timeInterval);

% SMA for TCS
Yk = zeros(length(Y),6);
for kk = 1:length(T)
    Yk(kk,1:6) = car2kep(Y(kk,1:3),Y(kk,4:6),data.mi);
end

vec = VectorLight-1;
vec = abs(vec);
idx = find(vec);
paraTCS = Yk(idx(1),:);

%% Sun vector for ADCS
% vecDates = data.InitDay +  T/86400;
% 
% vecSun = zeros(length(vecDates),3);
% for kk = 1:length(vecDates)
%    vecSun1 = uplanet(vecDates(kk),4); 
%    [vecSun(kk,:), ~] = kep2car(vecSun1,data.mi);
% end
% 
% incMars = 1.85061*pi/180;   % [rad]
% A = [cos(incMars) sin(incMars) 0
%      -sin(incMars) cos(incMars) 0
%      0 0 1];
% 
% Sun2Mars = zeros(length(T),3);
% for jj = 1:length(T)
%     Sun2Mars(jj,:) = A*vecSun(jj,1:3)';
% end
% 
% Mars2sc = Y(:,1:3);
% 
% SunDirection = - (Mars2sc + Sun2Mars);

%% Graphs
close all

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

plot3(Y(:, 1), Y(:, 2), Y(:, 3),'Color',[0.9490, 0.4745, 0.3137]), hold on, axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')

% Thrust profile
figure()
%sgtitle("Propulsion system on/off", 'FontSize', 20)
subplot(3,1,1)
plot(T/86400, abs(VectorThrust), 'Color',[0.1020, 0.6667, 0.74120])
ylabel('Thrust [N]','FontSize',20)
hold on
plot([T(1) T(end)]/86400, [0.1509 0.1509], 'r')
hold off
xlim([T(1) T(end)/86400])

% On/Off propulsion system
subplot(3,1,2)
plot(T/86400, VectorLight, 'Color',[155/255, 155/255, 155/255])
yticks([0 1])
yticklabels({'false','true'})
ylabel('Eclipse','FontSize',20)
xlim([T(1) T(end)/86400])

subplot(3,1,3)
plot(T/86400, VectorThRange, 'Color',[0.9490, 0.4745, 0.3137])
xlabel('T [days]')
yticks([0 1])
yticklabels({'false','true'})
ylabel('$\frac{de}{dt} > 0$','Interpreter','Latex','FontSize',20)
xlim([T(1) T(end)/86400])

% RISULTATI NS THRUST
% NS1 = -0.1509;
% NS2 = -0.1513;
% NS3 = -0.1796;

% id = find(VectorThrust);
% nonZ = VectorThrust(id);
% VALUE = max(nonZ);
% 
% % [0.9490, 0.4745, 0.3137]
% % [0.1020, 0.6667, 0.74120]
% % [155/255, 155/255, 155/255]
% 
% figure()
% subplot(2,3,1)
% plot(T,Yk(:,1))
% ylabel('SMA')
% 
% subplot(2,3,2)
% plot(T,Yk(:,2))
% ylabel('ecc')
% 
% subplot(2,3,3)
% plot(T,Yk(:,3)*180/pi)
% ylabel('incl')
% 
% subplot(2,3,4)
% plot(T,Yk(:,4)*180/pi)
% ylabel('RAAN')
% 
% subplot(2,3,5)
% plot(T,Yk(:,5)*180/pi)
% ylabel('omeghino')
% 
% subplot(2,3,6)
% plot(T,Yk(:,6)*180/pi)
% ylabel('theta')



