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
  
CALLED FUNCTIONS:

AUTHORS: D'andrea Biagio | Frulla Piergiorgio | Inno Adriano Filippo

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