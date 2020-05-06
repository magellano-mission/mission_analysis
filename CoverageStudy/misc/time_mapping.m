function time_map = time_mapping(curr, YYY, T, THETA, h_user, data)

% Setup
TT = curr.walker(1);             % total number of satellites
Nplanes = curr.walker(2);              % number of planes
n_sat = TT/Nplanes;

n_lon = length(data.lon);
n_lat = length(data.lat);

time_map = zeros(n_lon, n_lat, data.NT);
for jj = 1 : data.NT
    instant_states = YYY(:, :, jj);
    instant_map = zeros(n_lon, n_lat);           % how many sats can view a certain point in a certain istant of time
    for lo = 1 : n_lon
        lon = data.lon(lo);
        for la = 1 : n_lat
            lat = data.lat(la);
            % for each point, each satellite...
            for ii = 1 : Nplanes
                for kk = 1 : n_sat
                    %for each timestep, checking if point P is inside footprint
                    kkk = coverageNumber(lon, lat, T(jj), instant_states((ii-1)*n_sat + kk, :),...
                        THETA((ii-1)*n_sat + kk, jj), h_user((ii-1)*n_sat + kk, jj), data);
                    instant_map(lo, la) = kkk + instant_map(lo, la);
                end
            end
        end
    end
    time_map(:,:,jj) = instant_map;
    
end

end
