function  r = kep2car_r_only(OrbPar)
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

%}

% taking the input variables
a = OrbPar(1);                       
e = OrbPar(2);                   
i = OrbPar(3);                 
RAAN = OrbPar(4);                    
PA = OrbPar(5);               
theta = OrbPar(6);               

% computing the position vector
if e >= 1
    p = a*(e^2 - 1);
else
    p = a*(1 - e^2);
end
norm_r = p./(1+e*cos(theta));                 
r_pf = norm_r.*[cos(theta); sin(theta); 0];     
[R] = body2perif_rot(i, RAAN, PA);          
r = R' * r_pf;                              

end

function   [R] = body2perif_rot(i, RAAN, PA)

%{
body2perif_rot.m - This function returns the rotation matrix from the
                   centered inertial reference frame to the perifocal 
                   reference frame.

INPUTS:
    i  [1]  the Inclination of the orbit [rad]
    RAAN   [1]  the Right Ascension of the Ascending Node of the orbit [rad]
    PA     [1]  the Anomaly of the Pericentre of the orbit [rad]
  
OUTPUTS:
    R    [3x3]  the rotation matrix       

---------------------------------------------------------------------------------------------------
%}

% computing the three rotation matrices
R_RAAN = [  cos(RAAN)    sin(RAAN)     0  ;
           -sin(RAAN)    cos(RAAN)     0  ;
                0           0          1 ];
                                
R_INCLI =     [      1       0             0       ;
                     0   cos(i)    sin(i)  ;
                     0  -sin(i)    cos(i) ];

R_PA = [  cos(PA)   sin(PA)   0  ;
         -sin(PA)   cos(PA)   0  ;
            0         0       1] ;
                
R = R_PA*R_INCLI*R_RAAN;

end
