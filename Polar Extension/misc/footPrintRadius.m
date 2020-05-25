function [theta, h_sat, h_user] = footPrintRadius(gamma, Y, data)
% gamma: bw/2 [deg]
% Y: position state [km]
% altitude: [km] 

% Initialization
n = length(Y(:,1));
h_sat = zeros(n,1);
theta = zeros(n,1);
mask = data.mask*pi/180;

% Conversion
for kk = 1 : n                              % here probable error if alt ~0, because it won't compute satellite altitude
    if data.alt ~= 0
        [h_sat, h_user] = MCI2altitude(Y(kk,:));
        h_user = h_user + data.alt;
    else
        [h_sat, h_user] = MCI2altitude(Y(kk,:));
    end
    
    theta(kk) = - mask + abs(acos(cos(mask)/(1 + h_sat / h_user)));
    h_sat(kk) = h_sat;
    theta_max_inspection = rad2deg(theta(kk))
end

end