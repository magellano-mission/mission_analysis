function a_SRP = SRP(t, R_orbit, data)
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

CR = data.CR;                           % reflection coefficient
Acs = data.Acs;                         % lateral aerea
m_sat = data.m_sat;                     % satellite mass
InitData = data.initialDate;            % initial date
R_pl = data.R_pl;                       % planet radius
P0 = 4.5*1e-6;                          % solar radiation pressure coefficient
AU = 1.496*1e8;                         % astronomic unit
AOverM = Acs/m_sat;                     % Am parameter

InitDay = date2mjd2000(InitData);       % initial day
DayActual = InitDay + t/86400;

[r_S] = - uplanet(DayActual, 4);        % retriving the sun position vector wrt the planet
light = los(R_orbit, r_S, R_pl);        % checking the sight of ligth

if light
    rSun2Sc = R_orbit - r_S;
    a_srp = P0 * AU^2 / ((norm(rSun2Sc))^2) * CR * AOverM;
    a_SRP = a_srp * rSun2Sc / norm(rSun2Sc) / 1e3; % convert to km/s^2
else
    a_SRP = [0 0 0]';
end