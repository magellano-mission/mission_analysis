function a_Moon = MoonPerturb(DayActual, R_orbit, data, Moon)
%{
MoonPerturb.m - This function returns the perturbing acceleration due to
                the effects of the Moon as a third body.

INPUTS:
    DayActual   [1] the time in MJD2000; [s]
    R_orbit   [3x1] the position vector of the spacecraft in the geocentric
                    reference frame; [km]
  
OUTPUTS:
    a_Moon [3x1]  the vector of the perturbing acceleration in geocentric reference 
                  frame. [km/s^2]      
  
%}

mi_moon = data.moon;
if isrow(R)
    R_orbit = R';
end

switch Moon
    
    case 'Phobos'
        
        
    case 'Deimos'
        
        
    case 'Earth_moon'
    [R_moon, ~] = ephMoon(DayActual); % position of moon in inertial frame
    
end

if isrow(R_moon)
    R_moon = R_moon';
end

r_moon = norm(R_moon);
rSc2Moon = R_moon - R_orbit; % relative vector between satellite and moon
r_rel = norm(rSc2Moon);

% Vector of perturbing acceleration of moon%% Moon's perturbation
a_Moon = mi_moon*((rSc2Moon)./(r_rel^3) - R_moon./(r_moon^3)); % [km/s^2]

end