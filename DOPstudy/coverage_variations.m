% Coverage: Figures of Merit
close all; clear; clc

%% Adding to path
% addpath(genpath(fileparts(pwd)))

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
run('config.m')

%% check
check_errors(data);

%% Sim
tic
[Cov_Results, time_map, GDOP_map] = NestedFor(data);
toc

%% plot
sim_plots(data, Cov_Results, time_map);

T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
tspan = linspace(0, data.N_orbits*T_orb, data.NT);
figure()
plot(tspan/3600,squeeze(GDOP_map(10,10,:)))
