function [lat, lon, rad] = groundTrack(r, t)
%r(j) : column vector 3x1
%We have to include the case of keplerian parameters instead of cartesian
%ones.

    % Greenwich initial position
    thetaG = 0;     %[rad]
    % Earth rotation rate
    T = 1.026 * 24 * 3600;
    we = 2 * pi / T;    %[rad/s]
    
[alpha, delta] = angles(r);
lat = delta;
lon = alpha - (thetaG*ones(size(r,1),1) + we*t);

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