function  [r, v] = kep2car(OrbPar)
%{
geo2kep.m - This function returns the radius and the velocity in a certain
            position on the orbit knowing the orbital parameters of the orbit 
            and the gravitational constant of the celestial body.

INPUTS:
    obital_parameters    [6x1]-structure with this ordered variables:
                                  - a    [1]   semimajor axis of the orbit; [km]
                                  - e  [1]   eccentricity of the orbit;
                                  - i [1]  Inclination of the orbit; [deg]
                                  - RAAN  [1]  Right Ascension of the
                                               Ascending Node of the orbit; [deg]
                                  - PA    [1]  Argument of the Pericentre
                                               of the orbit; [deg] 
                                  - theta [1]  argument of the point on the
                                               orbit. [deg]
    mi       [1]  the gravitational constant of the celestial body in one
                  of the foci of the orbit. [km^3/s^2]
  
OUTPUTS:
    r_geo  [3x1]  the position of the orbiting body in the centered
                  inertial reference frame; [km]
    v_geo  [3x1]  the velocity of the orbiting body in the centered
                  inertial reference frame; [km/s]
    
CALLED FUNCTIONS:  -  body2perif_rot.m

AUTHORS: D'andrea Biagio | Frulla Piergiorgio | Inno Adriano Filippo

---------------------------------------------------------------------------------------------------
%}

% taking the input variables
a = OrbPar(1);                       
e = OrbPar(2);                   
i = OrbPar(3);                 
RAAN = OrbPar(4);                    
PA = OrbPar(5);               
theta = OrbPar(6);               

% computing the position vector
p = a*(1-e^2);                                
norm_r = p./(1+e*cos(theta));                 
r_pf = norm_r.*[cos(theta); sin(theta); 0];     
[R] = body2perif_rot(i, RAAN, PA);          
r = R' * r_pf;                              

end