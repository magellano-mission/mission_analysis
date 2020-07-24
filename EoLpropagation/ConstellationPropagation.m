% End of life propagation
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

%% post-process
Y = const.Y;
Ycar = zeros(length(T),3);
for k = 1:length(T)
    Ycar(k,:) = kep2car(Y(k,:),data.mi);
end

%% Plots
plot3(Ycar(:,1),Ycar(:,2),Y(:,3))
axis equal
