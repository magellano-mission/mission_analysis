% Coverage: Figures of Merit
close all; clear; clc

%% Adding to path
addpath(genpath(fileparts(pwd)))

%% Figure Initialization
set(0,'DefaultFigureUnits', 'normalized');
set(0,'DefaultFigurePosition',[0 0 1 1]);
set(0,'DefaultTextFontSize',18);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultAxesXGrid','on')
set(0,'DefaultAxesYGrid','on')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

%%
run('data.m')

%% check
check_errors(SimType, wlk_vec, inclinations, bw, semi_major_axes);

%% Sim
[Min_cov_lat, time_map, N_min, N_mean] = NestedFor(wlk_vec, inclinations, bw, semi_major_axes, lon, lat, timesteps, N_orbits, alt, para, perturb);
save('matriciozza.mat', 'Min_cov_lat')

%% plot
sim_plots(SimType, wlk_vec, inclinations, bw, semi_major_axes, Min_cov_lat, lon, lat, N_min);
