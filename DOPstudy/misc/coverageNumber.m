function [N, rho] = coverageNumber(lon, lat, t, Y, theta, h_user, data)

% Lon, lat of the target [deg]
% t 
% Y state vector
% theta: footprint area angle @ center [rad]
% h of the user

Y_unit = Y/norm(Y);
[r, posMCI] = lla2MCI(lon, lat, t, h_user, data);      % target position vector 

% Check if target central angle is lower than maximum Mars central angle
if acos(dot(r, Y_unit)) < (theta)
    N = 1;                           % inside the geometric horizon
    rho = processForQ(Y, posMCI);

else
    N = 0;                        % outside the geometric horizon
    rho = zeros(size(Y));
end

 
end



%{
Come sistemare:
- FATTO calcolare modulo di rho (vettore che va da O a B) in ECI (rimane uguale)
- FATTO girare la terna da MCI a TH fino a matchare con la TH centrata in O,
  quindi la rotazione che ho fatto dovrebbe essere giusta (con alpha e delta 
  dello user) --> il centro è ancora in C (figura libro Curtis)
- spostare il centro da C a O (però essendo ellissoide R non è allineato
  lungo z, perché z parte da C' --> capire come fare)
- calcolare A e a come faccio sulle mie dispense a pag. 50
- usare questi angoli per trovare B in coordinate TH
- ora la posizione dello user è 0,0,0
- proseguire come fatto fino ad ora per creare distM e Q
%}
