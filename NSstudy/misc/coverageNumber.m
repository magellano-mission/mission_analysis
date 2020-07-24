function [N, rho] = coverageNumber(lon, lat, t, Y, theta, h_user, data)

% Lon, lat of the target [deg]
% t 
% Y state vector
% theta: footprint area angle @ center [rad]
% h of the user

Y_unit = Y/norm(Y);
[r, posMCI] = lla2MCI(lon, lat, t, h_user, data);    % target position vector 

% Check if target central angle is lower than maximum Mars central angle
if acos(dot(r, Y_unit)) < (theta)
    N = 1;                           % inside the geometric horizon
    rho = Y - posMCI';   
else
    N = 0;                           % outside the geometric horizon
    rho = zeros(size(Y));
end

end