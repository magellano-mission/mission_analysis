% Main
close all; clear; clc

%% Adding to path
addpath(genpath(fileparts(strcat(pwd, '/misc'))))

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
run('config.m')

tic
[T, Ysat] = ode113(@GaussIntegration, [0, data.tend], data.X0_kep, data.opt, data);
t1 = toc

% Vector containing al the time instants in mjd2000
vectorDates = data.InitDay + T;
MarsIncl = 1.85061*pi/180;
rotMatrix = [1 0 0
             0 cos(MarsIncl) sin(MarsIncl)
             0 -sin(MarsIncl) cos(MarsIncl)];

ra = zeros(length(vectorDates),1);
dec = zeros(length(vectorDates),1);

tic
for kk = 1:length(vectorDates)
    MarsPosKep = uplanet(vectorDates(kk),4);
    MarsPosCart = kep2car(MarsPosKep, data.mi);
    MarsPosCart = - rotMatrix*MarsPosCart;
    [ra(kk), dec(kk)] = RaDec_from_r(MarsPosCart);
end
t2 = toc

beta = asin(cos(dec).*sin(Ysat(:,3)).*sin(Ysat(:,4)-ra) + sin(dec).*cos(Ysat(:,3)));

%%
plot(T,beta*180/pi)

