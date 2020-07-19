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


if isrow(R_orbit)
    R_orbit = R_orbit';
end

switch Moon
    
    case 'Phobos'
        mi_moon = data.miPhob;
        [Rmoon] = MoonEphs(DayActual, data, Moon);
        
    case 'Deimos'
        mi_moon = data.miDeim;
        [Rmoon] = MoonEphs(DayActual, data, Moon);
        
    case 'Earth_moon'
    [Rmoon, ~] = ephMoon(DayActual); % position of moon in inertial frame
    
end

if isrow(Rmoon)
    Rmoon = Rmoon';
end

r_moon = norm(Rmoon);
rSc2Moon = Rmoon - R_orbit; % relative vector between satellite and moon
r_rel = norm(rSc2Moon);

% Vector of perturbing acceleration of moon%% Moon's perturbation
a_Moon = mi_moon*((rSc2Moon)./(r_rel^3) - Rmoon./(r_moon^3)); % [km/s^2]

end