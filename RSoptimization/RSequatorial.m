clc; clear; close all;
load('beamwith.mat')

hmin = 500;
RM = 3393;
h_sat = SMa - RM;
SMa = SMa(h_sat > hmin);
Psat = P_sat((h_sat > hmin), :)';
h_sat = SMa - RM;

Nb = length(bw);
Na = length(SMa);
Ptot = zeros(Nb, Na);
Nsat = zeros(Nb, Na);
R = zeros(Nb, Na);
for i = 1:Nb
    gamma = deg2rad(bw(i)/2);
    for j = 1:Na
        R(i, j) = abs(pi/2 - gamma - acos( cos(pi/2 - gamma) * (1 + h_sat(j) / RM) ));
        Nsat(i, j) = ceil(pi/R(i, j)); 
        Ptot(i, j) = Psat(i, j)*Nsat(i, j);
    end
end

Rdeg = R*180/pi;
LatRange = 15;

feasible_map = Rdeg >= sqrt(2)*LatRange;


% inf_map = not(feasible_map);
Ptot_feasible = feasible_map.*Ptot;
Ptot_feasible(not(feasible_map)) = NaN;
[Pmin, rows] = min(Ptot_feasible);
[Pmin, col] = min(Pmin);
row = rows(col); 
amin = SMa(col);
bmin = bw(row);
Nsatmin = Nsat(row, col);
P_per_sat = Pmin/Nsatmin;


surf(bw, SMa, Ptot', 'EdgeColor', 'none')
xlabel('bw [deg]'), ylabel('sma [km]'), zlabel('Ptot')
title('equatorial RS')

clearvars -except amin bmin Nsatmin Pmin P_per_sat LatRange
