% Coverage: Figures of Merit
close all; clear; clc

%% Adding to path
%addpath(genpath(fileparts(pwd)))

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
[Cov_Results, time_map, GDOP_map, state] = NestedFor(data);
toc

%% plot
close all
sim_plots(data, Cov_Results, time_map);

% curiosity coord: lon = 137.4 lat = -4.6
figure()
subplot(3,1,1)
T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
tspan = linspace(0, data.N_orbits*T_orb, data.NT);
plot(tspan/3600, squeeze(GDOP_map(159,86,:)), 'Color', 'k')
ylabel('GDOP')
legend('Curiosiry','Location','northwest')
grid minor

% viking 1 coord: lon -50 lat 22.5
subplot(3,1,2)
T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
tspan = linspace(0, data.N_orbits*T_orb, data.NT);
plot(tspan/3600, squeeze(GDOP_map(65,113,:)), 'Color', 'k')
ylabel('GDOP')
legend('Viking 1','Location','northwest')
grid minor

% Skipper: lon = -167 lat = -81
subplot(3,1,3)
T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
tspan = linspace(0, data.N_orbits*T_orb, data.NT);
plot(tspan/3600, squeeze(GDOP_map(8,10,:)), 'Color', 'k')
xlabel('Time [hr]')
ylabel('GDOP')
legend('Skipper','Location','northwest')
grid minor

