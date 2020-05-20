function [value,isterminal,direction] = event_interp_LT(t_seconds, x, parameters)

%{ %event function that stops integration when s/c reaches SOI 
%  INPUT : 
% 
% % x : [6,1] planetorcentric states (CARTESIAN COORDS)
% % parameters : parameters of the s/c
% % t : seconds
% % % %  OUTPUT 
% pos_diff : distance from SOI
% isterminal : switch to halt integration 
% direction : direction of the states to reach the event
%}

    %parameteres :    
    isterminal = ones(2,1); %stops integration if 1
   
    muS = astroConstants(4);
    muP = astroConstants(14);
    
    %computation of conditions
    r_sc = x(1:3);
    r_sc_norm = norm(r_sc);
    kepPlanet = uplanet(t_seconds/86400,4);
    rP = kep2car2(kepPlanet, muS);
    rP_norm = norm(rP);
    r_scP = r_sc - rP; %distance s-c planet
    r_scPnorm = norm(r_scP);
    
    %SOI definition
    rSOI = rP_norm * (muP/muS)^(2/5);

    value = [r_sc_norm - rP_norm; ...
             rSOI - r_scPnorm]; 
         
    %isterminal defied previously
    direction = [0; 0];
end

