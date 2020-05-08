function [c, ceq] = mycon(x)

%%%%%%% constraints %%%%%%%%
 P_user_max = 30;
 LatRange = 25;
 
 %%%%%% data %%%%%%%%
 mods_user = 1;
 data_rate_user = 4e6;
 N = 1;  
 %f = 401e6;
 f = 8.5e9; 
 lambda = 3e8/f;
 RM = 3389.5; 
 T_eq = 250;         %[K]
 k = 1.3806e-23;     %[?]

% Semimajor axis:
SMA = x(1);            %[km]

% Reciever Gain:
G_user = x(2);            %[dBi]

% Satellite beamwidth:
bw = x(3);  %[deg]


G_sat = 10 * log10(30000 / bw^2);
bandwidth_sat = data_rate_user ./ mods_user * N; 

if mods_user == 1
    CNRs_user = 10.6;
elseif mods_user == 2
    CNRs_user = 13.6;
end
                               
h = SMA - RM;
losses = 20 * log10(4 * pi * h*1e3 / lambda);          %[dB]
design_margin = 3;                                  %[dB]

P_noise_sat = 10 * log10(k * T_eq .* bandwidth_sat);  

P_rec_sat = CNRs_user + P_noise_sat; 

P_user = P_rec_sat - G_sat - G_user + losses + design_margin; 

gamma = deg2rad(bw/2);
Rdeg = abs(pi/2 - gamma - acos( cos(pi/2 - gamma) * (1 + h / RM) ))*180/pi;

c(1) = P_user - P_user_max;
c(2) = sqrt(2)*LatRange - Rdeg;

ceq = 0;

end
