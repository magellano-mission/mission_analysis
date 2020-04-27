clear, clc, close all

mu = astroConstants(14);
R = 3389.5;

n = 2000;

SMA = 11500;
INC = deg2rad(55);

% SMA = 20700;
% INC = deg2rad(0);

X0 = [SMA, 1e-8, INC, 0, 0, 0];

T = 2 * pi * sqrt(SMA^3 / mu);
tspan = linspace(0, 6*T, n);

gamma = deg2rad(20);

opt = odeset('AbsTol', 1e-10, 'RelTol', 1e-10);

[T,Y] = trajectory(X0, tspan, mu, opt);

[lat, lon] = groundTrack(Y, T);
theta = footPrintRadius(gamma, Y, 0);

n = length(T);

%%
figure
plot(lon, lat, 'r.', 'LineWidth', 1)
hold on

for kk = 220 : 220
    plotCoverage(Y(kk, :), T(kk), theta(kk))
end

grid on
xline(0);
yline(30, '--');
yline(-30, '--');
yline(0);
xlim([-180 180])
ylim([-90 90])

%%
N = zeros(n);

for kk = 220 : 220
    N(kk) = coverageNumber(20, 50, T(kk), Y(kk,:), theta(kk));
end

%%
figure
plot(T, N, 'k', 'LineWidth', 1)
ylim([0 2])
grid on
