function [time_map, GDOP_map] = time_mapping(curr, YYY, T, THETA, h_user, data)

% Setup
TT = curr.walker(1);                   % total number of satellites
Nplanes = curr.walker(2);              % number of planes
n_sat = TT/Nplanes;                    % saellites per plane

if data.PE == true
    n_sat_pol = data.Nsat_pol;
end

n_lon = length(data.lon);
n_lat = length(data.lat);

time_map = zeros(n_lon, n_lat, data.NT);
GDOP_map = zeros(n_lon, n_lat, data.NT);

study = data.study;

switch study
    case "coverage"
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
                            % for each timestep, checking if point P is inside footprint
                            kkk = coverageNumber(lon, lat, T(jj), instant_states((ii-1)*n_sat + kk, :),...
                                THETA((ii-1)*n_sat + kk, jj), h_user((ii-1)*n_sat + kk, jj), data);
                            instant_map(lo, la) = kkk + instant_map(lo, la);
                        end
                    end
                    if data.PE == true
                    for pol = 1 : n_sat_pol
                         kkk = coverageNumber(lon, lat, T(jj), instant_states(Nplanes*n_sat + pol, :),...
                                THETA(Nplanes*n_sat + pol, jj), h_user(Nplanes*n_sat + pol, jj), data);
                         instant_map(lo, la) = kkk + instant_map(lo, la);
                    end
                    end
                end
            end
            time_map(:,:,jj) = instant_map;

        end
        
    case "DOP"
        for jj = 1 : data.NT
            distM = zeros(TT, 3);
            instant_states = YYY(:, :, jj);
            instant_map = zeros(n_lon, n_lat);           % how many sats can view a certain point in a certain istant of time
            GDOP_inst_map = zeros(n_lon, n_lat);

            for lo = 1 : n_lon
                lon = data.lon(lo);                     % [deg]
                for la = 1 : n_lat          
                    lat = data.lat(la);                 % [deg]
                    % for each point, each satellite...
                    for ii = 1 : Nplanes
                        for kk = 1 : n_sat
                            %for each timestep, checking if point P is inside footprint
                            [kkk, rho] = coverageNumber(lon, lat, T(jj), instant_states((ii-1)*n_sat + kk, :),...
                                THETA((ii-1)*n_sat + kk, jj), h_user((ii-1)*n_sat + kk, jj), data);
                            instant_map(lo, la) = kkk + instant_map(lo, la);
                            if sum(rho) ~= 0
                                distM(instant_map(lo, la),:) = rho/norm(rho);
                            end
                        end
                    end
                    if data.PE == true
                    for pol = 1 : n_sat_pol
                         [kkk, rho] = coverageNumber(lon, lat, T(jj), instant_states(Nplanes*n_sat + pol, :),...
                                THETA(Nplanes*n_sat + pol, jj), h_user(Nplanes*n_sat + pol, jj), data);
                         instant_map(lo, la) = kkk + instant_map(lo, la);
                            if sum(rho) ~= 0
                                distM(instant_map(lo, la),:) = rho/norm(rho);
                            end
                    end
                    end
                    distM2 = [distM, -ones(TT,1)];
                    distM2(all(distM2 == 0,2),:) = [];
                    Q = inv(distM2'*distM2);
                    GDOP_inst_map(lo,la) = sqrt(trace(Q));
                end
            end

            time_map(:,:,jj) = instant_map;
            GDOP_map(:,:,jj) = GDOP_inst_map;
            GDOP_map(time_map < 4) = 20;

        end

end
