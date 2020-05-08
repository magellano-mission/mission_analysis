% Station Keeping Box
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
run('ETconfig.m')

lb = 5e5;
ub = 2e6;

T = -1e-3;

%% Boundary refinition
% mi = astroConstants(14);
% r0 = 10500;
% m = 250;
% dth = pi/3;

% [lb, ub] = ETphasBoundaries(lb, ub, mi, r0, m, T, dth);

%% optimization
Ttot = EToptim(T, lb, ub, data); % total time of the phasing [days]

