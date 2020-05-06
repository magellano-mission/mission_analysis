function [lat, lon, rad] = groundTrack(r, t)
%r(j) : column vector 3x1
%We have to include the case of keplerian parameters instead of cartesian
%ones.

    % Airy0 crater initial position
    theta_Airy_0 = 0;     %[rad]
    % Mars rotation rate
    T = 1.026 * 24 * 3600;
    wM = 2 * pi / T;    %[rad/s]
    
[alpha, delta] = angles(r);
lat = delta;
lon = alpha - (theta_Airy_0*ones(size(r,1),1) + wM*t);

rad = [lat, lon];
while any(lon > pi)
    lon(lon>pi) = lon(lon>pi) - 2*pi;
end

while any(lon < -pi)
    lon(lon<-pi) = lon(lon<-pi) + 2*pi;
end
lat = rad2deg(lat);
lon = rad2deg(lon);
end