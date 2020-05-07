% Station Keeping Box
close all; clear; clc

%% Adding to path
addpath(genpath(fileparts(strcat(pwd, '/functions'))))

%% loading phobos ephemerides
load('PhobosEphs.mat')

%% Figure Initialization
set(0,'DefaultFigureUnits', 'normalized');
set(0,'DefaultFigurePosition', [0 0 1 1]);
set(0,'DefaultTextFontSize', 18);
set(0,'DefaultAxesFontSize', 18);
set(0,'DefaultAxesXGrid', 'on')
set(0,'DefaultAxesYGrid', 'on')
set(0,'defaultLegendInterpreter', 'latex');
set(0,'defaultAxesTickLabelInterpreter', 'latex');

%%
run('config_propagation.m')

%% sim
tic
const = PropagatedOrbits(data);
toc

%% SK maneuver
dv = SKmaneuver([const.Y(end, 1), const.Y(end, 2), const.Y(end, 3)], data);

%% plots

figure, plot(const.T, const.Y(:, 1)), ylabel('a');
figure, plot(const.T, const.Y(:, 2)), ylabel('e');
figure, plot(const.T, rad2deg(const.Y(:, 3))), ylabel('i');
figure, plot(const.T, rad2deg(const.Y(:, 4))), ylabel('RAAN');
figure, plot(const.T, rad2deg(const.Y(:, 5))), ylabel('PA');
figure, plot(const.T, rad2deg(const.Y(:, 6))), ylabel('theta');
% figure, plot3(const.Y(:, 1), const.Y(:, 2), const.Y(:, 3))

