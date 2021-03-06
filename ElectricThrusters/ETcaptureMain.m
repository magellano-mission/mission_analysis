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
TT = 5e7;                                % long period, Event function inside the ode
data.T = -0.2;                           % [N] for clearness the Thrust is here
[T, Y] = ode113(@ETcaptIntegration, [0, TT], data.Y0, data.opt1, data);
kep = car2kep(Y(end, 1:3), Y(end, 4:6), data.mi);
toc

%% Post-process
r_fin = norm(Y(end,1:3))

% K = zeros(length(T),6);
% 
% for k = 1:length(T)
%     K(k,:) = car2kep(Y(k,1:3),Y(k,4:6), data.mi);
% end

a_J2 = kep(1);
e_GNSS = kep(2);
i_GNSS = kep(3);
K = -3*data.J2*sqrt(data.mi)*data.R_pl^2./(2*a_J2.^(7/2)*(1-e_GNSS^2)^2);
RAAN_dot = K*cos(i_GNSS);       % variation per seconds in rad

rpM_orbit = 206644545;  % [km]
eM_orbit = 0.09341233;
smaM_orbit = rpM_orbit/(1 - eM_orbit);
wM_orbit = sqrt(astroConstants(4)/smaM_orbit^3);

% delta T between arrival of the first stack and the second to get in
% correct RAAN
deltaT = (120*pi/180)/(wM_orbit - RAAN_dot)/86400;

%% Graphs
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
plot3(Y(:, 1), Y(:, 2), Y(:, 3)), axis equal, hold on

hold off
figure()
plot(T,Y(:,7))
    
    
    
    
    
    
