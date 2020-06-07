function [T] = d2T (d,data)

% Function to compute maximum available thrust at a given distance from the Sun

% Input 
% d : distance from the Sun [km]

% Output
% T : maximum thrust available [N]

% Parameters
AU = 149597870.7;   % 1 AU [km]
A = data.Apanels;        % Solar panels surface [m2] for NS
% A = 42; % for ECS and RS

% Total available power
P = DistArray(d/AU, A);

% Available power for propulsion
P = P/1.2 - 220;

% Thrust from power
T = 3.042027115182780e-6*P.^2 + 0.022507117026406*P + 12.682081572453510; %[mN]
T = T/1000; % [N]

if T > 0.25
    T = 0.25;
end

%% Function to compute total available power
function[P] = DistArray(d, A)
    
% Input 
% d : distance from the Sun [km]
% A : solar panels surface [m2]

% Output
% P : maximum thrust available [N]

Xd = 0.85;          %Efficiency of the paths from the solar arrays directly to the loads 

SCe    = 0.316;     % solar cell efficiency [-]
P0     = irrSun(d); % Irradiance [W/m2]
Id     = 0.77;      % Inherent Degradation [-]
theta  = 5;         % Angle between Sun direction and normal to the panels [deg]
Degyr  = 0.5/100;   % Degradation per year [-]
LT     = 6;         % Lifetime of the mission [years]

Po   = P0*SCe;             % Output from solar cell type [W/m2]

PBol = Po*Id*cosd(theta);  % Beginning of Life Power [W/m2]

Ld   = (1- Degyr)^LT;      % Lifetime Degradation [-]

PEol = PBol *Ld;           % End of Life Power [W/m2]

Psa = A*PEol;              % Ideal solar Array Power [W]

P = Psa*Xd;                % Actual solar Array Power [W]
end

%% Function to compute solar irradiance
function[irr] = irrSun(d)
    
% Input 
% d : distance from the Sun [AU]

% Output
% irr : irradiance [W/m2]

c = 299792458;         % Speed of light [m/s]

P0 =  4.5e-6;          % Solar Radiation Pressure @ 1AU [Pa]

SRP = P0*(1/d)^2;      % Solar Radiation Pressure [Pa]

irr = SRP*c;           % Irradiance [W/m2]

end

end