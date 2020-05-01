function N = coverageNumber(lat, lon, t, Y, theta, h)
%lat, lon of the target
%t 
%Y state vector
%theta :footprint area angle @ center [deg]
Y = squeeze(Y);
Y = Y/norm(Y);
r = lla2MCI(lat, lon, t, h);

if acos(dot(r,Y)) < (theta)
    N = 1;
else
    N = 0;
end

 
end
