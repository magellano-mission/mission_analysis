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
run('ET_SKconfig.m')                       %%% ADD YOUR DATA

%%

TT = 1e8;                                  % long period, Event function inside the ode
data.T = 1e-1;                             % [N] for clearness the Thrust is here

[T, Y] = ode113(@ETIntegration, [0, TT], data.Y0, data.opt, data);
kep1 = car2kep(Y(end, 1:3), Y(end, 4:6), data.mi);
TOF = T(end)
plot3(Y(:, 1), Y(:, 2), Y(:, 3))

