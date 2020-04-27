function cov = getMinCoverage(walker, bw, SMA, INC, lon, lat, disc, alt)
% walker: [TT, P, F]
% bw: beamwidth (degree)
% SMA: semi-major axis [km]
% INC: inclination [degree]
% lon: [min_long max_long] (degree)
% lat: [min_lat max_lat] (degree)
% disc: [discret_lon discret_lat]
% alt: altitude (km)

if nargin == 7
    alt = 0;
end

mu = astroConstants(14);

TT = walker(1);
P = walker(2);
F = walker(3);

n = 2000;

n_sat = TT/P;
n_orbits = P;

n_lon = disc(1);
n_lat = disc(2);

LON = linspace(lon(1), lon(2), n_lon);
LAT = linspace(lat(1), lat(2), n_lat);


T_orb = 2 * pi * sqrt(SMA^3 / mu);
tspan = linspace(0, 3*T_orb, n);

X0 = [SMA, 1e-8, INC, 0, 0, 0];

gamma = deg2rad(bw/2);

opt = odeset('AbsTol', 1e-6, 'RelTol', 1e-6);

N = zeros(n,n_sat * n_orbits); %[time istant x number of satellite]
N_mesh = zeros(n_lat, n_lon);
N_mesh_mean = zeros(n_lat, n_lon);


for lo = 1 : n_lon
    for la = 1 : n_lat
        
        for ii = 1 : n_orbits
            for kk = 1 : n_sat
                X0(6) = (kk - 1) * pi /2 + (ii - 1) * F * 2 * pi / TT;
                X0(4) = (ii - 1) * 2 * pi / P;
                [T,Y] = trajectory(X0, tspan, mu, opt);
                theta = footPrintRadius(gamma, Y, alt);
                for jj = 1 : n
                    N(jj, (ii-1)*n_sat + kk) = coverageNumber(LAT(la), LON(lo), T(jj), Y(jj,:), theta(jj));
                end
            end
        end
        
        NN = sum(N,2);
        N_min = min(NN);
        N_mesh(la, lo) = N_min;
        N_mean = mean(NN);
        N_mesh_mean(la, lo) = N_mean;
    end
end

cov = min(N_min);

end