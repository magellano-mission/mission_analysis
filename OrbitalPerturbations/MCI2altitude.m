function[alt,R_real] = MCI2altitude(r)
% Compute the altitude starting from the position in ECI frame 
% Earth Geoid model is assumed 
%
% INPUT : 
% r [3x1] : Position Vector in ECI Frame [km]
% alt     : Altitude [km]
% R_real  : Earth radius at a given latitude [km]
%
% Other m-files required :  none 


r = r * 1000; % km-m conversion 
x = r(1); y = r(2) ; z = r(3); 

% WGS84 ellipsoid constants:
a = 3395428;
b = 3377678;
e = sqrt(1 - b^2 / a^2);

% calculations:
ep  = sqrt((a^2-b^2)/b^2);
p   = sqrt(x.^2+y.^2);
th  = atan2(a*z,b*p);

lat = atan2((z+ep^2.*b.*sin(th).^3),(p-e^2.*a.*cos(th).^3));
R_real   = a./sqrt(1-e^2.*sin(lat).^2);
alt = p./cos(lat)-R_real;

% Return Output in km 
alt = alt/1000; % km conversion 
R_real = R_real/1000;

return