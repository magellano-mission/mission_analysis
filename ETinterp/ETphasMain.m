% ET phasing 
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
run('ETphasConfig.m')                       %%% ADD YOUR DATA

%% Boundary refinition
lb = 5e5;                                   % define the boundaries
ub = 2e6;

T = -1e-3;                                  % [N] for clearness the Thrust is here

% mi = astroConstants(14);
% r0 = 10500;
% m = 250;
% dth = pi/3;

% [lb, ub] = ETphasBoundaries(lb, ub, mi, r0, m, T, dth);

%% optimization
Ttot = ETphasOptim(T, lb, ub, data);        % total time of the phasing [days]

