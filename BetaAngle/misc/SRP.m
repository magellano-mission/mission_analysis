function a_srp = SRP(t, R_orbit, data)
%{
SRP.m - This function returns the perturbing acceleration due to
        the effects of the Solar Radiation Pressure.

INPUTS:
    t         [1]  the time in MJD2000; [s]
    R_orbit [3x1]  the position vector of the spacecraft in the geocentric
                 reference frame; [km]
    data [17x1]-structure with these ordered variables:
                    - CR           [1]   the solar radiation pressure
                                         coefficient;
                    - Acs          [1]   the cross section of the spacecraft; [m^2]
                    - initialDate  [1x6] the date of the beginning of the
                                         simulation; [yyyy/mm/dd/hh/mm/ss]
                    - R_pl         [1]   the mean radius of the Earth; [km]
                    - mi           [1]   the gravitational constant of the
                                         Earth; [km^3/s^2]
                    - a0           [1]   semimajor axis of the given orbit; [km]
                    - e0           [1]   the eccentricity of the given
                                         orbit;
                    - i0           [1]   the inclination of the given orbit; [deg]
                    - m_sat        [1]   the mass of the spacecraft; [kg]
                    - J2           [1]   the value of the J2 perturbation
                                         parameter;
                    - beta         [1]   the ballistic coefficient;
                    - T            [1]   the period of the given orbit; [s]
                    - T_mult       [1]   the period of the repeating
                                         ground-track orbit; [s]
                    - a_mult       [1]   the semimajor axis of the
                                         repeating ground-track orbit; [s]
                    - tmax_unp     [1]   the last time of simulation for the 
                                         unperturbed orbit; [s]
                    - tmax_pert    [1]   the last time of simulation for the 
                                         perturbed orbit; [s]
                    - option       [1]   options for the ODE equation.
                  
  
OUTPUTS:
    a_SRP  [3x1]  the vector of the perturbing acceleration in geocentric reference 
                  frame. [km/s^2]      
  
CALLED FUNCTIONS: - date2mjd2000.m
                  - sunPosition.m
                  - los.m

AUTHORS: D'andrea Biagio | Frulla Piergiorgio | Inno Adriano Filippo

---------------------------------------------------------------------------------------------------
%}

CR = data.CR;                               % reflection coefficient
R_pl = data.R_pl;                           % planet radius
S0 = data.S0;                               % solar radiation pressure coefficient
eps = data.eps;                             
AOverM = data.AOverM;                       % Am parameter

DayActual = data.InitDay + t/86400;

KepS = uplanet(DayActual, 4);
rM2S = - kep2car_r_only(KepS);              % retriving the sun position vector wrt the planet
light = los(R_orbit, rM2S, R_pl);            % checking the sight of ligth

c = cos(eps); s = sin(eps);
E = [1 0 0
     0 c -s
     0 s c ];

rSun2Sc = R_orbit - rM2S;
r = rSun2Sc / norm(rSun2Sc);
S = S0/((norm(rSun2Sc))^2);
 
if light
    a_srp = S * CR * AOverM * r / 1e3;              % convert to km/s^2
else
    a_srp = [0 0 0]';
end

a_srp = E*a_srp;