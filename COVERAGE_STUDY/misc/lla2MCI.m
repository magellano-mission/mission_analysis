function [r] = lla2MCI(lat, lon, t, h, para)
 
lat = lat/180*pi;
lon = lon/180*pi;

curv_radius = para.rM_eq/sqrt(1-(sin(lat)*sin(para.ang_ecc))^2);  % curvature radius in the first vertical

% Conversion
delta = lat;
alpha = lon + (para.theta_Airy_0 + para.wM*t);
x = (curv_radius + h)*cos(delta)*cos(alpha);
y = (curv_radius + h)*cos(delta)*sin(alpha);
z = ((cos(para.ang_ecc)^2)*curv_radius + h)*sin(delta);
r = [x; y; z];
r = r/norm(r);

end

