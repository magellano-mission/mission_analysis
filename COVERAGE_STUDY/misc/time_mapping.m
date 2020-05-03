function time_map = time_mapping(walker, YYY, T, THETA, H, LON, LAT, para)

% Setup
timesteps = length(T);
TT = walker(1);             % total number of satellites
P = walker(2);              % number of planes
n_sat = TT/P;
n_orbits = P;
n_lon = length(LON);
n_lat = length(LAT);

time_map = zeros(n_lon, n_lat, timesteps);
for jj = 1 : timesteps
    instant_states = squeeze(YYY(:,jj,:));               % 2D matrix (dimension squeezed later)
    instant_map = zeros(n_lon,n_lat);           % how many sats can view a certain point in a certain istant of time
    for lo = 1 : n_lon
        for la = 1 : n_lat
            % for each point, each satellite...
            for ii = 1 : n_orbits
                for kk = 1 : n_sat
                    %for each timestep, checking if point P is inside footprint
                    kkk = coverageNumber(LAT(la), LON(lo), T(jj), instant_states((ii-1)*n_sat + kk,:), THETA((ii-1)*n_sat + kk,jj), H((ii-1)*n_sat + kk,jj), para);
                    instant_map(lo,la) = kkk + instant_map(lo,la);
                end
            end
        end
    end
    time_map(:,:,jj) = instant_map;
    
end

end
