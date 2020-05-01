function N = coverageNumber(lat, lon, t, Y, theta, h)
%lat, lon of the target
%t 
%Y state vector
%theta: footprint area angle @ center [deg]

Y = squeeze(Y);                   % satellite position vector at instant t
Y = Y/norm(Y);
r = lla2MCI(lat, lon, t, h);      % target position vector 

% Check if target central angle is lower than maximum Mars central angle
if acos(dot(r,Y)) < (theta)
    N = 1;                        % inside the geometric horizon
else
    N = 0;                        % outside the geometric horizon
end
 
end
