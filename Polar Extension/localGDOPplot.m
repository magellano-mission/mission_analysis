%% Load data
clc;clear;close all

load('NominalConditions')
GDOP_nominal = GDOP_map;

load('OffNominalConditions')
GDOP_offNominal = GDOP_map;

%% plot
close all

% curiosity coord: lon = 137.4 lat = -4.6
figure()
subplot(3,1,1)
T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
tspan = linspace(0, data.N_orbits*T_orb, data.NT);
plot(tspan/3600, squeeze(GDOP_nominal(159,86,:)), 'Color', [0.1020, 0.6667, 0.74120],'Linewidth',1.5)
hold on
plot(tspan/3600, squeeze(GDOP_offNominal(159,86,:)), 'Color', [0.1020, 0.6667, 0.74120],'LineStyle','-.','Linewidth',1.5)
hold off
ylabel('GDOP')
legend('Curiosity','Location','northwest')
grid minor


% viking 1 coord: lon -50 lat 22.5
subplot(3,1,2)
T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
tspan = linspace(0, data.N_orbits*T_orb, data.NT);
plot(tspan/3600, squeeze(GDOP_nominal(65,113,:)), 'Color', [155/255, 155/255, 155/255],'Linewidth',1.5)
hold on
plot(tspan/3600, squeeze(GDOP_offNominal(65,113,:)), 'Color', [155/255, 155/255, 155/255],'LineStyle','-.','Linewidth',1.5)
hold off
ylabel('GDOP')
legend('Viking 1','Location','northwest')
grid minor


% Skipper: lon = -167 lat = -81
subplot(3,1,3)
T_orb = 2 * pi * sqrt(data.sma^3 / data.mi);
tspan = linspace(0, data.N_orbits*T_orb, data.NT);
plot(tspan/3600, squeeze(GDOP_nominal(8,10,:)), 'Color', [0.9490,0.4745,0.3137],'Linewidth',1.5)
hold on
plot(tspan/3600, squeeze(GDOP_offNominal(8,10,:)), 'Color', [0.9490,0.4745,0.3137],'LineStyle','-.','Linewidth',1.5)
hold off
xlabel('Time [hr]')
ylabel('GDOP')
legend('Skipper','Location','northwest')
grid minor
