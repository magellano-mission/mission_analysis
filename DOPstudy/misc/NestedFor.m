function [Cov_Results, time_map, GDOP_map, state] = NestedFor(data)

for w = 1:data.Nw
    curr.walker = [data.Nsat(w), data.Norb, data.walk_phas];    % Walker vector [sats, planes, phase shift]
    
    for i = 1:data.Ni
        curr.inc = data.inc(i);                  % inclination [rad]
        
        for b = 1:data.Nb
            curr.beam = data.bw(b);              % beamwidth [deg]
            
                for s = 1:data.Na
                    curr.sma = data.sma(s);      % semi-major axis [km]
                    [YYY, T, THETA, ~, h_user] = const_orbits(curr, data);
                    [time_map, GDOP_map, state] = time_mapping(curr, YYY, T, THETA, h_user, data);
                    GDOP_max = max(GDOP_map, [], 3);
                    GDOP_max(GDOP_max == 100) = NaN;
%                     HDOP_max = max(HDOP_map, [], 3);
%                     HDOP_max(HDOP_max == 100) = NaN;
%                     VDOP_max = max(VDOP_map, [], 3);
%                     VDOP_max(VDOP_max == 100) = NaN;
                end
        end
    end
end

Cov_Results.GDOP_max = GDOP_max;
% Cov_Results.HDOP_max = HDOP_max;
% Cov_Results.VDOP_max = VDOP_max;

end
