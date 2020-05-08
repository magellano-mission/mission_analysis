function light_switch = los(Rsat, Rsun, R_pl)
% This function uses the ECI position vectors of the satellite (r_sat)
% and the sun (r_sun) to determine whether the earth is in the line of
% sight between the two.

%{
los.m - This function determines whether the spacecraft is seen
                   by the Sun or instead it is shadowed by the Earth.

INPUTS:
    Rsat   [3x1]  the position vector of the satellite wrt the planet; [km]
    Rsun   [3x1]  the position vector of the Sun wrt the planet; [km]
    R_pl   [1]    radius of the planet; [km]
  
OUTPUTS:
    light:switch  [1]  the switcher between "on light" position (true) and 
                       "in shadow" position (false);     
  

---------------------------------------------------------------------------------------------------
%}

rsat = norm(Rsat);
rsun = norm(Rsun);

% Angle between sun and satellite position vectors:
theta = acosd(dot(Rsat, Rsun)/rsat/rsun);

% Angle between the satellite position vector and the radial to the point
% of tangency with the earth of a line from the satellite:
theta_sat = acosd(R_pl/rsat);

% Angle between the sun position vector and the radial to the point
% of tangency with the earth of a line from the sun:
theta_sun = acosd(R_pl/rsun);

% Determine whether a line from the sun to the satellite
% intersects the earth:

if theta_sat + theta_sun <= theta
    light_switch = true; 
else
    light_switch = false; %no
end

end 