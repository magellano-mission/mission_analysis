function [r, posMCI] = lla2MCI(lon, lat, t, h_user, data)
 
lat = lat/180*pi;       % [rad]
lon = lon/180*pi;       % [rad]

curv_radius = data.rM_eq/sqrt(1-(sin(lat)*sin(data.ang_ecc))^2);  % curvature radius in the first vertical

% Conversion
delta = lat;            % [rad]
alpha = lon + (data.theta_Airy_0 + data.wM*t);                    % [rad]
x = (curv_radius + h_user)*cos(delta)*cos(alpha);
y = (curv_radius + h_user)*cos(delta)*sin(alpha);
z = ((cos(data.ang_ecc)^2)*curv_radius + h_user)*sin(delta);
posMCI = [x; y; z];
r = posMCI/norm(posMCI);

end

