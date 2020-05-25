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
sma = curr.sma;

% Orbits computation
YYY_wlk = zeros(TT, 3, data.NT);      % states matrix (3D)
THETA_wlk = zeros(TT, data.NT);        % matrix containing Earth central angles
h_sat_wlk = zeros(TT, data.NT);            % surface altitude of the geodetic point
h_user_wlk = zeros(TT, data.NT);  

 for ii = 1 : Nplane
     for kk = 1 : n_sat
         X0(6) = (kk - 1) * pi * 2 /n_sat + (ii - 1) * F * 2 * pi / TT;
         X0(4) = (ii - 1) * 2 * pi / Nplane;
         [T, Y] = trajectory(X0, tspan, data, opt, sma);
         [theta, h_sat_wlk((ii-1)*n_sat + kk, :), h_user_wlk((ii-1)*n_sat + kk, :)] = footPrintRadius(gamma, Y, data);
         YYY_wlk((ii-1)*n_sat + kk, :, :) = Y';
         THETA_wlk((ii-1)*n_sat + kk, :) = theta;
     end
 end

if data.PE == true
    % Polar extension
    YYY_pol = zeros(data.Nsat_pol, 3, data.NT);      % states matrix (3D)
    THETA_pol = zeros(data.Nsat_pol, data.NT);       % matrix containing Earth central angles
    h_sat_pol = zeros(data.Nsat_pol, data.NT);       % surface altitude of the geodetic point
    h_user_pol = zeros(data.Nsat_pol, data.NT);  

    X02 = [data.sma_pol, 1e-8, data.inc_pol, 0, 0, 0];       % initial condition

    n_sat_pol = data.Nsat_pol;
    sma_pol = data.sma_pol;
    for kk = 1 : n_sat_pol
            X02(6) = (kk - 1) * pi * 2 /n_sat_pol;
            [T, Y] = trajectory(X02, tspan, data, opt, sma_pol);
            [theta, h_sat_pol(kk, :), h_user_pol(kk, :)] = footPrintRadius(gamma, Y, data);
            YYY_pol(kk, :, :) = Y';
            THETA_pol(kk, :) = theta;
    end

    YYY = [YYY_wlk; YYY_pol];
    THETA = [THETA_wlk; THETA_pol];
    h_sat = [h_sat_wlk; h_sat_pol];
    h_user = [h_user_wlk; h_user_pol];

else

    YYY = YYY_wlk;
    THETA = THETA_wlk;
    h_sat = h_sat_wlk;
    h_user = h_user_wlk;
end

    

end