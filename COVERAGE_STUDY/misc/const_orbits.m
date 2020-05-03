function [YYY, T, THETA, h_sat, h_user] = const_orbits(walker, bw, SMA, INC, timesteps, N_orbits, alt, perturb)

if nargin == 6
    alt = 0;
end

mu = astroConstants(14);            % Mars planetary constant

TT = walker(1);                     % total number of satellites
P = walker(2);                      % number of planes
F = walker(3);                      % relative spacing between satellites in adjacent planes
    

n_sat = TT/P;                       % satellites per plane
n_orbits = P;


T_orb = 2 * pi * sqrt(SMA^3 / mu);
tspan = linspace(0, N_orbits*T_orb, timesteps);

% Setup 
X0 = [SMA, 1e-8, INC, 0, 0, 0];     % initial condition
gamma = deg2rad(bw/2);              % satellite FoV
opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-13);

% Orbits computation
YYY = zeros(n_orbits* n_sat, timesteps, 3);      % states matrix (3D)
THETA = zeros(n_sat*n_orbits, timesteps);        % matrix containing Earth central angles
h_sat = zeros(n_sat*n_orbits, timesteps);            % surface altitude of the geodetic point
h_user = h_sat;

 for ii = 1 : n_orbits
     for kk = 1 : n_sat
         X0(6) = (kk - 1) * pi * 2 /n_sat + (ii - 1) * F * 2 * pi / TT;
         X0(4) = (ii - 1) * 2 * pi / P;
         [T, Y] = trajectory(X0, tspan, mu, opt, SMA, perturb);
         [theta, h_sat((ii-1)*n_sat + kk, :), h_user((ii-1)*n_sat + kk, :)] = footPrintRadius(gamma, Y, alt);
         YYY((ii-1)*n_sat + kk, :, :) = Y;
         THETA((ii-1)*n_sat + kk, :) = theta;
     end
 end

end
