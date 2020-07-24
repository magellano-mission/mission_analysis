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

%% delta quantities
T = const.T;
da = const.Y(:, 1) - data.X0_kep(1);
de = const.Y(:, 2) - data.X0_kep(2);
di = rad2deg(const.Y(:, 3) - data.X0_kep(3));
dO = rad2deg(const.Y(:, 4) - data.X0_kep(4));
do = rad2deg(const.Y(:, 5) - data.X0_kep(5));
dth = rad2deg(const.Y(:, 6) - data.X0_kep(6));

%% retriving worst points for the SK
[maxda, ia] = max(abs(da));
[maxde, ie] = max(abs(de));
[maxdi, ii] = max(abs(di));

%% SK maneuver
dv = SKmaneuver([const.Y(ia, 1), const.Y(ie, 2), const.Y(ii, 3)], data);  %[m/s]

%% report subplot

figure;
subplot(3, 1, 1), plot(T/86400, da) , ylabel('\Deltaa [km]'), xlabel('time [days]');
subplot(3, 1, 2), plot(T/86400, de), ylabel('e'), xlabel('time [days]');
subplot(3, 1, 3), plot(T/86400, di), ylabel('\Deltai [deg]'), xlabel('time [days]');
