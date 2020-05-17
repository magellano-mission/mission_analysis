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
[Cov_Results, Maps_Results] = NestedFor(data);
toc

switch data.study
    case "coverage"
        time_map = Maps_Results.time_map;
    case "DOP"
        GDOP_map = Maps_Results.GDOP_map;
        time_map = Maps_Results.time_map;
end

%% plot
close all
sim_plots(data, Cov_Results, time_map);

% curiosity coord: lon = 137.4 lat = -4.6
figure()
subplot(2,1,1)
T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
tspan = linspace(0, data.N_orbits*T_orb, data.NT);
plot(tspan/3600, squeeze(GDOP_map(88,34,:)), 'Color', 'k')
%xlabel('Time [hr]')
ylabel('GDOP')
legend('Curiosiry','Location','northwest')
grid minor
ylim([1 2.6])

% viking 1 coord: lon -50 lat 22.5
subplot(2,1,2)
T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
tspan = linspace(0, data.N_orbits*T_orb, data.NT);
plot(tspan/3600, squeeze(GDOP_map(37,44,:)), 'Color', 'k')
xlabel('Time [hr]')
ylabel('GDOP')
legend('Viking 1','Location','northwest')
grid minor
ylim([1 2.6])

% Skipper: lon = -167 lat = -81
figure()
T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
tspan = linspace(0, data.N_orbits*T_orb, data.NT);
plot(tspan/3600, squeeze(GDOP_map(5,4,:)), 'Color', 'k')
xlabel('Time [hr]')
ylabel('GDOP')
legend('Skipper','Location','northwest')
grid minor

% Colori Magellano
% [0.1020    0.6667    0.74120]
% [0.9490    0.4745    0.3137]