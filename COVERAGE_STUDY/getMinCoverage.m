function [N_min, N_mean, cov] = getMinCoverage(disc, time_map)
% walker: [TT, P, F]
% bw: beamwidth (degree)
% SMA: semi-major axis [km]
% INC: inclination [degree]
% lon: [min_long max_long] (degree)
% lat: [min_lat max_lat] (degree)
% disc: [discret_lon discret_lat]
% alt: altitude (km)

n_lon = disc(1);
n_lat = disc(2);

% N = zeros(timesteps,n_sat * n_orbits); %[time istant x number of satellite]
N_min = zeros(n_lon, n_lat);
N_mean = zeros(n_lon, n_lat);

    for lo = 1 : n_lon
        for la = 1 : n_lat
            N_min(lo,la) = min(time_map(lo,la,:));% N_min = min(NN); % N_mesh(la, lo) = N_min;
        end
    end
     for lo = 1 : n_lon
        for la = 1 : n_lat
            N_mean(lo,la) = mean(time_map(lo,la,:));% N_min = min(NN); % N_mesh(la, lo) = N_min;
        end
    end
N_mean = mean(time_map,3);% N_mean = mean(NN); %  N_mesh_mean(la, lo) = N_mean;
cov = min(N_min);
end