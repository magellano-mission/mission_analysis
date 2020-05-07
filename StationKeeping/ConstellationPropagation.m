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
dv = SKmaneuver([const.a(end), const.e(end), const.i(end)], data);

%% plots

figure, plot(const.T, const.a), ylabel('a');
figure, plot(const.T, const.e), ylabel('e');
figure, plot(const.T, const.i), ylabel('i');
figure, plot(const.T, const.RAAN), ylabel('RAAN');
figure, plot(const.T, const.PA), ylabel('PA');
figure, plot(const.T, const.th), ylabel('theta');
figure, plot3(const.Y(:, 1), const.Y(:, 2), const.Y(:, 3))

