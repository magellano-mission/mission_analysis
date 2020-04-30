function [r] = lla2MCI(lat, lon, t)

lat = deg2rad(lat);
lon = deg2rad(lon);

T = 1.02749125 * 24 * 3600;
wM = 2 * pi / T;    %[rad/s]
thetaG = 0;

delta = lat;
alpha = lon + (thetaG + wM*t);
x = cos(delta)*cos(alpha);
y = cos(delta)*sin(alpha);
z = sin(delta);
r = [x; y; z];
end

