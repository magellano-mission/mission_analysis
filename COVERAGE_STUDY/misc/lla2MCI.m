function [r] = lla2MCI(lat, lon, t, h)

rM_eq = 3393.4;         % equatorial radius [km]
rM_pol = 3375.7;        % polar radius [km]

ang_ecc = 2*atan(sqrt((rM_eq-rM_pol)/(rM_eq+rM_pol)));

lat = deg2rad(lat);
lon = deg2rad(lon);

curv_radius = rM_eq/sqrt(1-(sin(lat)*sin(ang_ecc))^2);

T = 1.02749125 * 24 * 3600;
wM = 2 * pi / T;    % [rad/s]
theta_Airy_0 = 0;   % Mars principal meridian

delta = lat;
alpha = lon + (theta_Airy_0 + wM*t);
x = (curv_radius + h)*cos(delta)*cos(alpha);
y = (curv_radius + h)*cos(delta)*sin(alpha);
z = ((cos(ang_ecc)^2)*curv_radius + h)*sin(delta);
r = [x; y; z];

end

