% ET SK 
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
run('ETcaptureConfig.m')                       %%% ADD YOUR DATA

%%

TT = 2e7;                                  % long period, Event function inside the ode
data.T = -2;                              % [N] for clearness the Thrust is here
[T1, Y1] = ode113(@ETcaptIntegration, [0, TT], data.Y0, data.opt1, data);
kep1 = car2kep(Y1(end, 1:3), Y1(end, 4:6), data.mi);

% TT = 3e6;                                % long period, Event function inside the ode
% data.T = 0;                              % [N] for clearness the Thrust is here
% [T2, Y2] = ode113(@ETIntegration, [0, TT], Y1(end, :), data.opt1, data);
% kep1 = car2kep(Y1(end, 1:3), Y1(end, 4:6), data.mi);

% data.T = -3;                              % [N] for clearness the Thrust is here
% [T2, Y2] = ode113(@EThypIntegration, [0, TT], Y1(end, :), data.opt2, data);
% kep2 = car2kep(Y2(end, 1:3), Y2(end, 4:6), data.mi);

% data.T = -1.2;                              % [N] for clearness the Thrust is here
% [T3, Y3] = ode113(@ETIntegration, [0, TT], Y2(end, :), data.opt3, data);
% kep3 = car2kep(Y3(end, 1:3), Y3(end, 4:6), data.mi);

plot3(Y1(:, 1), Y1(:, 2), Y1(:, 3)), axis equal, hold on
plot3(0, 0, 0, 'bo', 'MarkerSize', 15)
% plot3(Y2(:, 1), Y2(:, 2), Y2(:, 3))
% plot3(Y3(:, 1), Y3(:, 2), Y3(:, 3))

    
    
    
    
    
    
