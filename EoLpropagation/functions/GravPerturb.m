function [a_J2, a_J3, a_J4] = GravPerturb(t, R, data)
%{
GravPerturb.m - This function returns the perturbing acceleration due to
              the effects of Planet oblatness.

INPUTS:
    t       [1]  the time passed since the beginning of the simulation; [s]
    R     [3x1]  the position vector of the spacecraft in the geocentric
                 reference frame; [km]
 
                  
%}

R_pl = data.R_pl;                          % planet radius
J2 = data.J2;                              % J2 perturbation constant
J3 = data.J3;                              % J2 perturbation constant
J4 = data.J4;                              % J2 perturbation constant
mi = data.mi;                              % earth gravitational constant
w = data.wM;                                
r = norm(R);


% Rotation matrix from ECI to ECEF frame
B = [cos(w*t) sin(w*t) 0 ;
    -sin(w*t) cos(w*t) 0 ;
        0         0      1];

% S/c position in ECEF frame
s_sc = B*R;

% J2 Effect
k = 3/2 * J2 * mi * R_pl^2 / r^4;

a_j2_x = k*(s_sc(1)/r *(5 *s_sc(3)^2/r^2 -1));
a_j2_y = k*(s_sc(2)/r *(5 *s_sc(3)^2/r^2 -1));
a_j2_z = k*(s_sc(3)/r *(5 *s_sc(3)^2/r^2 -3));

a_J2 = [a_j2_x; a_j2_y; a_j2_z];                         % [km/s^2]

% J3 Effect
kz = 1/2 * J3 * mi * R_pl^3 / r^5;

a_j3_x = kz*5*s_sc(1)/r*(7*(s_sc(3)/r)^3-3*(s_sc(3)/r));
a_j3_y = kz*5*s_sc(2)/r*(7*(s_sc(3)/r)^3-3*(s_sc(3)/r));
a_j3_z = kz*3*(1-10*(s_sc(3)/r)^2+35/3*(s_sc(3)/r)^4);

a_J3 = [a_j3_x; a_j3_y; a_j3_z];                        % [km/s^2]

% J4 Effect
kj = 5/8 * J4 * mi * R_pl^4 / r^6;
a_j4_x = kj*(3 - 42*(s_sc(3)/r)^2 + 63*(s_sc(3)/r)^4)*s_sc(1)/r;
a_j4_y = kj*(3 - 42*(s_sc(3)/r)^2 + 63*(s_sc(3)/r)^4)*s_sc(2)/r;
a_j4_z = kj*(15 - 70*(s_sc(3)/r)^2 + 63*(s_sc(3)/r)^4)*s_sc(1)/r;

a_J4 = [a_j4_x; a_j4_y; a_j4_z];                        % [km/s^2]

% rotation in ECI
a_J2 = B'*a_J2;
a_J3 = B'*a_J3;
a_J4 = B'*a_J4;


end