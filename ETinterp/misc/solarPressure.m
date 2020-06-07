function a_srp = solarPressure(rr_sc, M, parameters)
% Function used to compute the perturbing acceleration due to the Solar
% radiation pressure. Line of sight and eclipse computation is included.

% INPUT   :  
% rr_sc   : [3x1] [km]  Spacecraft position in ECI frame
% data    : struct 
% t       : integration time 


% OUTPUT  : 
% a_srp   : [3x1] [km/s^2] Perturbing Acceleration 

c_r = parameters.c_r;
A_m = parameters.Across_sun / M;

AU = 1.496e8;

R = norm(rr_sc);

p_sr = 4.5e-6;

Asrp = p_sr * (AU^2) * c_r * A_m / R^2 * 1e-3;

a_srp = - Asrp / R * rr_sc; % [km / s^2]

end