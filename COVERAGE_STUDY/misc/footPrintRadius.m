function [theta, h] = footPrintRadius(gamma, Y, altitude)
%gamma: bw/2 [deg]
%Y : position state [km]
%altitude: [km] 
n = length(Y(:,1));

h = zeros(n,1);
theta = zeros(n,1);

for kk = 1 : n
    
    [h, R] = MCI2altitude(Y(kk,:));
    if nargin == 3
        R = R + altitude;
    end
    %theta = lambda0
    %gamma = eta
    theta(kk) = abs(pi/2 - gamma - acos( cos(pi/2 - gamma) * (1 + h / R) ));
    
    h(kk) = h;
end

end