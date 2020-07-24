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

%% Retrieve data
switch data.study
    case "coverage"
        time_map = Maps_Results.time_map;
    case "DOP"
        GDOP_map = Maps_Results.GDOP_map;
        time_map = Maps_Results.time_map;
end

%% Save data
if data.Nsat == 21 && data.study == "DOP"
    save('NominalConditions')
elseif data.Nsat == 18 && data.study == "DOP"
    save('OffNominalConditions')
end

%% plot
close all
sim_plots(data, Cov_Results, time_map);

if data.study == "coverage"
    figure()
    subplot(1,2,1)
    timeGap = Cov_Results.max_cov_gap;
    maximum = max(timeGap);
    maximum = maximum - maximum(end);
    gap = barh(data.lat, maximum);
    gap.FaceColor = [0.9490,0.4745,0.3137];
%     gap.EdgeColor = [1.949    0.4745    0.3137];
    ylabel('Latitude')
    xlabel('Maximum time gap [hr]')

    subplot(1,2,2)
    percCov = Cov_Results.percent_cov;
    minimum = min(percCov);
    perc = barh(data.lat, minimum);
    perc.FaceColor = [0.9490,0.4745,0.3137];
%     perc.EdgeColor = [0.9490    0.4745    0.3137];
    xlim([75 100])
    ylabel('Latitude')
    xlabel('Percent coverage [%]')
end

if data.study == "DOP"
    % curiosity coord: lon = 137.4 lat = -4.6
    figure()
    subplot(3,1,1)
    T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
    tspan = linspace(0, data.N_orbits*T_orb, data.NT);
    plot(tspan/3600, squeeze(GDOP_map(159,86,:)), 'Color', [0.1020, 0.6667, 0.74120])
    ylabel('GDOP')
    legend('Curiosiry','Location','northwest')
    grid minor
  

    % viking 1 coord: lon -50 lat 22.5
    subplot(3,1,2)
    T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
    tspan = linspace(0, data.N_orbits*T_orb, data.NT);
    plot(tspan/3600, squeeze(GDOP_map(65,113,:)), 'Color', [155/255, 155/255, 155/255])
    ylabel('GDOP')
    legend('Viking 1','Location','northwest')
    grid minor
   

    % Skipper: lon = -167 lat = -81
    subplot(3,1,3)
    T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
    tspan = linspace(0, data.N_orbits*T_orb, data.NT);
    plot(tspan/3600, squeeze(GDOP_map(8,10,:)), 'Color', [0.9490,0.4745,0.3137])
    xlabel('Time [hr]')
    ylabel('GDOP')
    legend('Skipper','Location','northwest')
    grid minor
end

