function [YYY, T, THETA, H] = const_orbits(walker, bw, SMA, INC, alt, timesteps, N_orbits)

if nargin == 7
    alt = 0;
end

mu = astroConstants(14);

TT = walker(1);
P = walker(2);
F = walker(3);


n_sat = TT/P;
n_orbits = P;


T_orb = 2 * pi * sqrt(SMA^3 / mu);
tspan = linspace(0, N_orbits*T_orb, timesteps);

X0 = [SMA, 1e-8, INC, 0, 0, 0];

gamma = deg2rad(bw/2);

opt = odeset('AbsTol', 1e-6, 'RelTol', 1e-6);


%orbits computation
YYY = zeros(n_orbits* n_sat, timesteps, 3);
THETA = zeros(n_sat*n_orbits, timesteps);
H = zeros(n_sat*n_orbits, timesteps);
 for ii = 1 : n_orbits
            for kk = 1 : n_sat
%                 X0(6) = (kk - 1) * pi /2 + (ii - 1) * F * 2 * pi / TT;
                X0(6) = (kk - 1) * pi * 2 /n_sat + (ii - 1) * F * 2 * pi / TT;
                X0(4) = (ii - 1) * 2 * pi / P;
                [T,Y] = trajectory(X0, tspan, mu, opt);
                [theta, h] = footPrintRadius(gamma, Y, alt);
                YYY((ii-1)*n_sat + kk,:,:) = Y;
                THETA((ii-1)*n_sat + kk,:) = theta;
                H((ii-1)*n_sat + kk,:) = h;
            end
 end

end
