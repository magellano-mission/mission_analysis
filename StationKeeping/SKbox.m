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

%%
run('config_box.m')

%% Sim
tic
[Cov_Results, time_map] = DummyPropagation(data);
toc

[flag] = CheckCoverage(Cov_Results, data);

if not(flag)
    [t_box] = ReconstructTime(time_map, data);
    [box] = ReconstructBox(t_box, data);
end



