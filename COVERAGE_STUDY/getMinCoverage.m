function [N_min, N_mean, cov] = getMinCoverage(disc, time_map)

% disc: [discret_lon discret_lat]
% time_map: temporal-spatil history

n_lon = disc(1);
n_lat = disc(2);

% Initialization
N_min = zeros(n_lon, n_lat);
N_mean = zeros(n_lon, n_lat);

% Computation of minimum
for lo = 1 : n_lon
    for la = 1 : n_lat
        N_min(lo,la) = min(time_map(lo,la,:));% N_min = min(NN); % N_mesh(la, lo) = N_min;
    end
end

% Computation of mean
for lo = 1 : n_lon
    for la = 1 : n_lat
        N_mean(lo,la) = mean(time_map(lo,la,:));% N_min = min(NN); % N_mesh(la, lo) = N_min;
    end
end

N_mean = mean(time_map,3);
cov = min(N_min);

end