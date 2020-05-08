function N = coverageNumber(lon, lat, t, Y, theta, h, data)
%lat, lon of the target
%t 
%Y state vector
%theta: footprint area angle @ center [deg]

Y_unit = Y/norm(Y);
r = lla2MCI(lon, lat, t, h, data);      % target position vector 

% Check if target central angle is lower than maximum Mars central angle
if acos(dot(r, Y_unit)) < (theta)
    N = 1;                        % inside the geometric horizon
else
    N = 0;                        % outside the geometric horizon
end
 
end
