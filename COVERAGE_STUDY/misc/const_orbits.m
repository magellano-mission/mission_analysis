function [YYY, T, THETA, h_sat, h_user] = const_orbits(curr, data)

TT = curr.walker(1);                     % total number of satellites
Nplane = curr.walker(2);                      % number of planes
F = curr.walker(3);                      % relative spacing between satellites in adjacent planes

n_sat = TT/Nplane;                       % satellites per plane

T_orb = 2 * pi * sqrt(curr.sma^3 / data.mi);
tspan = linspace(0, data.N_orbits*T_orb, data.NT);

% Setup 
X0 = [curr.sma, 1e-8, curr.inc, 0, 0, 0];       % initial condition
gamma = deg2rad(curr.beam/2);                   % satellite FoV
opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-13);

% Orbits computation
YYY = zeros(TT, data.NT, 3);      % states matrix (3D)
THETA = zeros(TT, data.NT);        % matrix containing Earth central angles
h_sat = zeros(TT, data.NT);            % surface altitude of the geodetic point
h_user = zeros(TT, data.NT);  

 for ii = 1 : Nplane
     for kk = 1 : n_sat
         X0(6) = (kk - 1) * pi * 2 /n_sat + (ii - 1) * F * 2 * pi / TT;
         X0(4) = (ii - 1) * 2 * pi / Nplane;
         [T, Y] = trajectory(X0, tspan, data, opt, curr);
         [theta, h_sat((ii-1)*n_sat + kk, :), h_user((ii-1)*n_sat + kk, :)] = footPrintRadius(gamma, Y, data);
         YYY((ii-1)*n_sat + kk, :, :) = Y;
         THETA((ii-1)*n_sat + kk, :) = theta;
     end
 end

end
