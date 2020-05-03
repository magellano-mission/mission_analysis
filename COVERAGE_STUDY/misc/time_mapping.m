<<<<<<< HEAD
function time_map = time_mapping(walker, YYY, T, THETA, H, LON, LAT, para)
=======
function [time_map, LON, LAT] = time_mapping(walker, YYY, T, THETA, H, lon, lat, disc)
>>>>>>> parent of e2edd09... abbiamo fatto casino (ma in meglio)

% Setup
timesteps = length(T);
TT = walker(1);             % total number of satellites
P = walker(2);              % number of planes
n_sat = TT/P;
n_orbits = P;
<<<<<<< HEAD
n_lon = length(LON);
n_lat = length(LAT);
=======

n_lon = disc(1);            % number of discretizing points in longitude
n_lat = disc(2);            % number of discretizing points in latitude

% delta_lon = (lon(2) - lon(1))/n_lon;    % longitude step
% delta_lat = (lat(2) - lat(1))/n_lat;    % latitude step
% 
% % Discertization
% LON = (lon(1)+delta_lon/2:delta_lon:lon(2)-delta_lon/2);
% LAT = (lat(1)+delta_lat/2:delta_lat:lat(2)-delta_lat/2);

LON = linspace(lon(1), lon(2), n_lon);
LAT = linspace(lat(1), lat(2), n_lat);
>>>>>>> parent of e2edd09... abbiamo fatto casino (ma in meglio)

time_map = zeros(n_lon, n_lat, timesteps);
for jj = 1 : timesteps
    instant_states = YYY(:,jj,:);               % 2D matrix (dimension squeezed later)
    instant_map = zeros(n_lon,n_lat);           % how many sats can view a certain point in a certain istant of time
    for lo = 1 : n_lon
        for la = 1 : n_lat
            % for each point, each satellite...
            for ii = 1 : n_orbits
                for kk = 1 : n_sat
                    %for each timestep, checking if point P is inside footprint
                    kkk = coverageNumber(LAT(la), LON(lo), T(jj), instant_states((ii-1)*n_sat + kk,1,:), THETA((ii-1)*n_sat + kk,jj), H((ii-1)*n_sat + kk,jj));
                    instant_map(lo,la) = kkk + instant_map(lo,la);
                end
            end
        end
    end
    time_map(:,:,jj) = instant_map;
    
end

end
