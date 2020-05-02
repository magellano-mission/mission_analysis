% Coverage: Figures of Merit
close all; clear; clc
% load('matriciozza.mat')


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

%% Sim
tic
[Min_cov_lat, time_map, N_min, N_mean] = NestedFor(wlk_vec, inclinations, bw, semi_major_axes, lon, lat, timesteps, N_orbits, alt, para);
save('matriciozza.mat', 'Min_cov_lat')
toc
%% plot
sim_plots(SimType, wlk_vec, inclinations, bw, semi_major_axes, Min_cov_lat);
