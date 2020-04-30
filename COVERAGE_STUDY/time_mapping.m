function [time_map, LON, LAT] = time_mapping(walker, YYY, T, THETA, lon, lat, disc)

timesteps = length(T);

TT = walker(1);
P = walker(2);
F = walker(3);
n_sat = TT/P;
n_orbits = P;

n_lon = disc(1);
n_lat = disc(2);

delta_lon = (lon(2) - lon(1))/n_lon;
delta_lat = (lat(2) - lat(1))/n_lat;

%discertization
LON = (lon(1) + delta_lon/2:delta_lon:lon(2)-delta_lon/2);
LAT = (lat(1) + delta_lat/2:delta_lat:lat(2)-delta_lat/2);

time_map = zeros(n_lon, n_lat, timesteps);
for jj = 1 : timesteps
   waitforit = waitbar(jj/timesteps);
    instant_states = YYY(:,jj,:);  %2D matrix (dimension squeezed later)
    instant_map = zeros(n_lon,n_lat); %quanti satelliti vedono ciascun punto in un dato istante
    for lo = 1 : n_lon
        for la = 1 : n_lat
            %for each point, each satellite...
                   for ii = 1 : n_orbits
                   for kk = 1 : n_sat
                    %for each timestep, checking if point P is inside
                    %footprint
                    kkk = coverageNumber(LAT(la), LON(lo), T(jj), instant_states((ii-1)*n_sat + kk,1,:), THETA((ii-1)*n_sat + kk,jj));
                    instant_map(lo,la) = kkk + instant_map(lo,la);
                   end
                   end
        end
    end
time_map(:,:,jj) = instant_map; %  NN = sum(N,2);
%         %operation for each (lat,lon)
end
delete(waitforit)
end