
clc, clear, close all
h = 5000;

mu_m = astroConstants(14);
R_m = 3389.5;

a = R_m + h;
e = 0.00001;
i = deg2rad(50);
Om = 0;
om = 0;
theta = 0;

x_0 = [a, e, i, Om, om, theta];

T = 2 * pi * sqrt(a^3 / mu_m);
t_span = [0 100*T];

[t, x] = ode113(@gauss, t_span,  x_0);

m = length(t);

parout = zeros(m,3);

for kk = 1 : m
    [~, parout(kk,:)] = gauss(t(kk), x(kk,:));
end

parout_norm = parout./parout(:,1);

figure
for kk = 1 : 5
subplot(2,3,kk)
plot(t, x(:,kk))
grid on
hold on
plot(t, x(1,kk) * ones(length(t),1));
end

subplot(2,3,6)
semilogy(t, parout_norm);
grid on

figure
plot(x(:,3), x(:,4))
grid on

%%
state = x_0;
h_range = [200 6000];
n = 10;
n_orb = 100;

[kep_ranges] = ranges(state, h_range, n, n_orb);

INCs = rad2deg(kep_ranges(:,1:2));
RAANs = kep_ranges(:,3:4);

figure
for kk = 1 : n
xline(INCs(kk, 1), 'Color', kk / n *[1, 0, 0])
hold on
xline(INCs(kk, 2), 'Color', kk / n *[1, 0, 0])
yline(RAANs(kk, 1), 'Color', kk / n *[0, 1, 0])
yline(RAANs(kk, 2), 'Color', kk / n *[0, 1, 0])
end
grid on
xlim([min(min(INCs))*0.9999, max(max(INCs))*1.0001])
ylim([min(min(RAANs))*0.9999, max(max(RAANs))*1.0001])