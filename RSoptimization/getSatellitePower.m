function [Ptot] = getSatellitePower(x)
% P:                power needed for satellite transmitter [W]
%
%   x:  [1x3]       [SMA (km), G_user (dB), bw (deg)]
%                   
%   P_user_max:     maximum power for user transmitter [W]
%   f:              frequency needed [Hz]
%   mods_sat:       satellite downlink modulation 1 for BPSK, 2 for QPSK) []
%   mods_user:      user uplink modulation (1 for BPSK, 2 for QPSK) []
%   data_rate_sat:  satellite downlink data-rate [bps]
%   data_rate_user: user uplink data-rate [bps]
%   N:              number of users [ ]

 %f = 401e6;             % UHF
 f = 8.5e9;              % X-band 
 mods_sat = 1;
 data_rate_sat = 1e6;

% Semimajor axis:
SMA = x(1);            %[km]

% Reciever Gain:
G_user = x(2);            %[dBi]

% Satellite beamwidth:
bw = x(3);  %[deg]

% Wavelength definition
c = 3e8;            %[m/s]
lambda = c / f;     %[m]

% Path loss
RM = 3389.5;                                  %[m]
h = SMA - RM;
losses = 20 * log10(4 * pi * h*1e3 / lambda);          %[dB]
design_margin = 3;                                  %[dB]

% Noise power data
k = 1.3806e-23;     %[?]
T_eq = 250;         %[K]

% Data rate to achieve to the user:
bandwidth_user = data_rate_sat ./ mods_sat;      %[Hz]

% Minimum CNR to achieve:
if mods_sat == 1
    CNRs_sat = 10.6;
elseif mods_sat == 2
    CNRs_sat = 13.6;
end

% Noise Power:
P_noise = 10 * log10(k * T_eq .* bandwidth_user);        %[dBW]

% Minimum Recieved Power:
P_rec = CNRs_sat + P_noise;     %[dBW]

% Satellite EIRP Needed:
EIRP = P_rec - G_user + losses + design_margin;  %[dBW]

% Gain of the satellite computation wrt beamwidth:
G_sat = 10 * log10(30000 / bw^2);

% Power needed for the satellite:
Psat = 10.^((EIRP - G_sat)./10);

gamma = deg2rad(bw/2);              % half beamwidth
R = abs(pi/2 - gamma - acos( cos(pi/2 - gamma) * (1 + h / RM) ));   % radius [rad]
Nsat = ceil(pi/R);                  % minimum number of satellites
Ptot = Psat*Nsat;                   % total power

end
