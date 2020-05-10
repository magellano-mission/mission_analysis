function rho_vec = processForQ(Y, posMCI)
% DA CORREGGERE PER INSERIRE IL LOCAL HOUR ANGLE (DEFINIZIONE DI SIDEREAL
% TIME DI MARTE)

% lat = lat/180*pi;   % [rad]
% theta = data.theta_Airy_0 + data.wM*t + lon/180*pi;
% 
% Position vector of the satellite relative to the site
rho_vec = Y - posMCI';     
% rho = norm(rho_vec);
% 
% % Right ascension and declination of the satellite relative to the site (topocentric)
% [alpha, delta] = angles(rho_vec);
% 
% % Local hour angle
% t = theta - alpha;  % sbagliato, perché il 90 va da West a Sigma 
% 
% % Elevation and Azimuth
% h = asin(sin(lat)*sin(delta) + cos(lat)*cos(delta)*cos(t));
% A = - asin(sin(t)/cos(h)*cos(delta));
% 
% % Topocentric horizon
% rhoTH = rho*[cos(h)*sin(A), cos(h)*cos(A), sin(h)];



end