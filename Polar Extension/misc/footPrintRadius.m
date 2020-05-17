function [theta, h_sat, h_user] = footPrintRadius(gamma, Y, data)
% gamma: bw/2 [deg]
% Y: position state [km]
% altitude: [km] 

% Initialization
n = length(Y(:,1));
h_sat = zeros(n,1);
theta = zeros(n,1);

% Conversion
for kk = 1 : n                              % here probable error if alt ~0, because it won't compute satellite altitude
    if data.alt ~= 0
        h_user = data.R_eq + data.alt;
    else
        [h_sat, h_user] = MCI2altitude(Y(kk,:));
    end
    
    theta(kk) = abs(pi/2 - gamma - acos( cos(pi/2 - gamma) * (1 + h_sat / h_user) ));
    h_sat(kk) = h_sat;
end

end