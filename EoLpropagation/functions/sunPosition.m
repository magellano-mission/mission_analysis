function [r_S] = sunPosition(initialdate, AU, t)
%{
sunPosition.m - This function returns the position of the Sun with respect
                to the Earth.

INPUTS:
    initialdate  [1]  the date of the beginning of the simulation in
                      MJD2000; [s]
    AU           [1]  the Astronomic Unit; [km]
    t            [1]  the time since the beginning of the simulation; [s]

OUTPUTS:
    r_s        [3x1]  the position vector of the Sun in the geocentric reference frame; [km]      
  
CALLED FUNCTIONS:

AUTHORS: D'andrea Biagio | Frulla Piergiorgio | Inno Adriano Filippo

---------------------------------------------------------------------------------------------------
%}

% Initial day in MJD2000
n = initialdate;

% We add to initial day the days (fraction of days) from 
% the beginning of integration
n = n + t/(24*60*60); 

% Mean anomaly [deg]
M = 357.528 + 0.9856003*n; 
M = mod(M,360); 

% Mean longitude [deg]:
L = 280.460 + 0.98564736*n;
L = mod(L,360);

% Apparent ecliptic longitude [deg]  
lambda = L + 1.915*sind(M) + 0.020*sind(2*M); % from - Astronomical Almanac
lambda = mod(lambda,360);

% Obliquity of the ecliptic (deg):
eps = 23.439 - 0.0000004*n;

% Unit vector from earth to sun:
u = [cosd(lambda); sind(lambda)*cosd(eps); sind(lambda)*sind(eps)];

% Distance from earth to sun (km):
rS = (1.00014 - 0.01671*cosd(M) - 0.000140*cosd(2*M))*AU;

% Geocentric position vector (km):
r_S = rS*u;

end

