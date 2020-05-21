%Conway-based approach

%definition of parameters
close all, clear, clc
% Figure Initialization    
load('MagellanoColorMap.mat');
DefaultOrderColor = get(0, 'DefaultAxesColorOrder');
NewOrderColor = [0.9490    0.4745    0.3137
                 0.1020    0.6667    0.74120
                 155/255   155/255   155/255
                 DefaultOrderColor];  
             
set(0,'DefaultFigureColormap', MagellanoColorMap);
set(0, 'DefaultAxesColorOrder', NewOrderColor);
set(0,'DefaultLineLineWidth', 2)
set(0,'DefaultLineMarkerSize', 10)
set(0, 'DefaultFigureUnits', 'normalized');
set(0, 'DefaultFigurePosition', [0 0 1 1]);
set(0, 'DefaultTextFontSize', 18);
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

%initialization
data_stacks.t0sym = date2mjd2000([2025 1 1 0 0 0]);
data_stacks.tmax = date2mjd2000([2030 1 1 0 0 0]);

data_stacks.isPerturbed = 0;
data_stacks.isInterp = 1;
data_stacks.Isp = 4300;                                      % specific impulse [s]
data_stacks.Mdry = 8000;                                      % Total Mass of the s/c [kg]

[kepEarth, muS] = uplanet(data_stacks.t0sym, 3);

%definition of interplanetary arc
TOF = 3.5*365;%[days]
N_rev = 2; %max 3-4 for Conway hp

[r, TH, L, gamma1, gamma2, gamma, a_inplane, vr, vt, v1tra, v2tra, vnorm, time, dmdt, m, T_inplane, TOFr] = ...
Conway(TOF, N_rev, data_stacks);

%capture 
%%%%%%

%plots
if ~isnan(r) 
RE = zeros(length(TOFr), 3); RM = RE;
REnorm = zeros(length(TOFr),1); RMnorm = REnorm;
for i =1:length(TOFr)
    kepE = uplanet(data_stacks.t0sym+TOFr(i), 3);
    kepM = uplanet(data_stacks.t0sym+TOFr(i), 4);
    RE(i,:) = kep2car2(kepE, muS);
    RM(i,:) = kep2car2(kepM, muS);
    REnorm(i) = norm(RE(i,:));
    RMnorm(i) = norm(RM(i,:));
end

figure()
sgtitle('Sun distance')
plot(TOFr,RMnorm, 'DisplayName', 'Mars'), hold on
plot(TOFr,REnorm, 'DisplayName', 'Earth'), hold on
plot(TOFr, r, 'DisplayName', 's/c'), hold on, legend()

figure()
sgtitle('Thrust Profile')
subplot(3,1,1), plot(TOFr, a_inplane), title('a_{inplane}')
subplot(3,1,2), plot(TOFr, T_inplane), title('inplane T')
subplot(3,1,3), plot(TOFr, m), title('mass')

figure()
plot(r.*cos(TH), r.*sin(TH)), hold on, 
plot(norm(RE(1,:)),0,'o'), hold on
plot(norm(RM(end,:))*cos(L), norm(RM(end,:))*sin(L),'o'),axis equal

figure()
sgtitle('Flight path angle')
plot(TOFr, gamma), hold on 
yline(gamma1); hold on, yline(gamma2);

figure()
plot(TOFr, vt,'DisplayName','$v_{t}$'), hold on
yline(norm(v1tra),'DisplayName','vEt'); hold on, yline(norm(v2tra),'DisplayName','vMt');
plot(TOFr, vnorm,'DisplayName','$nomr(v)$'), hold on, plot(TOFr, vr, 'DisplayName','$v_{r}$')
legend()
else
    fprintf('No real solution for Conway algorithm \n')
end

