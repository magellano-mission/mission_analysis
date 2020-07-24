function [Cov_Results, Maps_Results] = NestedFor(data)

study = data.study;

switch study
    case "coverage"
        Min_cov_latlon = zeros(data.Na, data.Ni, data.Nw);
        Min_cov_lat = zeros(data.Na, length(data.lat));

        for w = 1:data.Nw
            curr.walker = [data.Nsat(w), data.Norb, data.walk_phas];

            for i = 1:data.Ni
                curr.inc = data.inc(i);             % inclination [rad]

                for s = 1:data.Na
                    curr.sma = data.sma(s);              % semi-major axis [km]
                    
                    [YYY, T, THETA, ~, h_user] = const_orbits(curr, data);
                    time_map = time_mapping(curr, YYY, T, THETA, h_user, data);
                    [Cov_Results] = getMinCoverage(time_map, data, T);
                    Min_cov_latlon(s, i, w) = Cov_Results.cov_mars;
                    Min_cov_lat(s, :) = Cov_Results.cov_lon';
                end
            end
        end
        
        Maps_Results.time_map = time_map;
        Cov_Results.Min_cov_latlon = Min_cov_latlon;
        Cov_Results.Min_cov_lat = Min_cov_lat;
        
    case "DOP"
        for w = 1:data.Nw
            curr.walker = [data.Nsat(w), data.Norb, data.walk_phas];    % Walker vector [sats, planes, phase shift]

            for i = 1:data.Ni
                curr.inc = data.inc(i);                  % inclination [rad]

                for s = 1:data.Na
                    curr.sma = data.sma(s);      % semi-major axis [km]

                    [YYY, T, THETA, ~, h_user] = const_orbits(curr, data);
                    [time_map, GDOP_map] = time_mapping(curr, YYY, T, THETA, h_user, data);
                    GDOP_max = max(GDOP_map, [], 3);
                    GDOP_max(GDOP_max == 20) = NaN;
                end
            end
        end

        Maps_Results.GDOP_map = GDOP_map;
        Maps_Results.time_map = time_map;
        Cov_Results.GDOP_max = GDOP_max;

end

end
