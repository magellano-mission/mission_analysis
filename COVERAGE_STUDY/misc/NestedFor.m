function [Cov_Results, time_map] = NestedFor(data)

Min_cov_latlon = zeros(data.Nb, data.Na, data.Ni, data.Nw);
Min_cov_lat = zeros(data.Nb, data.Na, length(data.lat));

for w = 1:data.Nw
    curr.walker = [data.Nsat(w), data.Norb, data.walk_phas];
    
    for i = 1:data.Ni
        curr.inc = data.inc(i);
        
        for b = 1:data.Nb
            curr.beam = data.bw(b);       % inclination [deg]
            
                for s = 1:data.Na
                    curr.sma = data.sma(s);              % semi-major axis [km]
                    
                    [YYY, T, THETA, ~, h_user] = const_orbits(curr, data);
                    time_map = time_mapping(curr, YYY, T, THETA, h_user, data);
                    [Cov_Results] = getMinCoverage(time_map, data, T);
                    Min_cov_latlon(b, s, i, w) = Cov_Results.cov_mars;
                    Min_cov_lat(b, s, :) = Cov_Results.cov_lon';
                end
        end
    end
end

Cov_Results.Min_cov_latlon = Min_cov_latlon;
Cov_Results.Min_cov_lat = Min_cov_lat;