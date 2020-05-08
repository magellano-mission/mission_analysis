function [r, posMCI] = lla2MCI(lon, lat, t, h, data)
 
lat = lat/180*pi;
lon = lon/180*pi;

curv_radius = data.rM_eq/sqrt(1-(sin(lat)*sin(data.ang_ecc))^2);  % curvature radius in the first vertical

% Conversion
delta = lat;
alpha = lon + (data.theta_Airy_0 + data.wM*t);
x = (curv_radius + h)*cos(delta)*cos(alpha);
y = (curv_radius + h)*cos(delta)*sin(alpha);
z = ((cos(data.ang_ecc)^2)*curv_radius + h)*sin(delta);
posMCI = [x; y; z];
r = posMCI/norm(posMCI);

end

