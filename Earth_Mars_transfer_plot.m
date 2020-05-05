function [e, Dv] = Earth_Mars_transfer_plot(Earth_time, Mars_time)

[k_E, mu_s] = uplanet(0, 3);
a_E = k_E(1);
T_E = 2*pi*sqrt(a_E^3/mu_s);
t_E = linspace(0, T_E, 500);
rr_E = zeros(3,length(t_E));

for i = 1:length(t_E)
    
    tE = t_E(i);
    [k_E, ~] = uplanet(tE/24/3600, 3);
    [r_E, ~] = kep2car2([k_E(1), k_E(2), k_E(3), k_E(4), k_E(5), k_E(6)], mu_s);
    rr_E(:, i) = r_E;
    
end

[k_M, ~] = uplanet(0, 4);
a_M = k_M(1);
T_M = 2*pi*sqrt(a_M^3/mu_s);
t_M = linspace(0, T_M, 500);
rr_M = zeros(3,length(t_M));

for i = 1:length(t_M)
    
    tM = t_M(i);
    [k_M, ~] = uplanet(tM/24/3600, 4);    
    [r_M, ~] = kep2car2([k_M(1), k_M(2), k_M(3), k_M(4), k_M(5), k_M(6)], mu_s);    
    rr_M(:, i) = r_M;
    
end

% COMPUTING TRANSFER ORBITS
%--------------------------------------------------------------------------

[k_E, ~] = uplanet(Earth_time, 3);
[r_E, v_E] = kep2car2([k_E(1), k_E(2), k_E(3), k_E(4), k_E(5), k_E(6)], mu_s);
        
[k_M, ~] = uplanet(Mars_time, 4);
[r_M, v_M] = kep2car2([k_M(1), k_M(2), k_M(3), k_M(4), k_M(5), k_M(6)], mu_s);
           
[~, ~, e, ~, VI, VF, ~, ~] = lambertMR(r_E, r_M, (Mars_time - Earth_time)*24*3600 , mu_s);

Dv = norm(VI' - v_E) + norm(v_M - VF');

y01 = [r_E; VI'];

[~, y1] = ode113( @(t, y) elliptical_orbit(t, y, mu_s) , linspace(Earth_time, Mars_time, 500)*24*3600, y01, odeset('RelTol', 1e-13, 'AbsTol', 1e-14));

E_d = mjd20002date(Earth_time);
E_date = strcat(num2str(E_d(3)), '/', num2str(E_d(2)), '/', num2str(E_d(1)));

M_d = mjd20002date(Mars_time);
M_date = strcat(num2str(M_d(3)), '/', num2str(M_d(2)), '/', num2str(M_d(1)));

% PLOTTING
%--------------------------------------------------------------------------

plot(rr_E(1, :), rr_E(2, :), '-', 'DisplayName', 'Earth orbit')
title('Interplanetary Transfer Orbit')
hold on
grid on

plot(rr_M(1, :), rr_M(2, :), '-', 'DisplayName', 'Mars orbit')

plot(y1(:, 1), y1(:, 2), 'Linewidth', 2, 'DisplayName', 'Transfer arc', 'Color', 'r')

text(r_E(1), r_E(2), E_date, 'Interpreter', 'Latex', 'FontSize', 8)
text(r_M(1), r_M(2), M_date, 'Interpreter', 'Latex', 'FontSize', 8)


legend('Location', 'South')
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 11)

ax = gca;
AU = 149597870.700;
ax.XAxis.TickValues = [-32:8:32]*AU;
ax.XAxis.TickLabels = {'-32', '-24', '-16', '-8', '0', '8', '16', '24', '32'};
ax.XLabel.String = 'x direction [AU]';
ax.XLabel.Interpreter = 'Latex';

ax.YAxis.TickValues = [-32:8:32]*AU;
ax.YAxis.TickLabels = {'-32', '-24', '-16', '-8', '0', '8', '16', '24', '32'};
ax.YLabel.String = 'y direction [AU]';
ax.YLabel.Interpreter = 'Latex';



end