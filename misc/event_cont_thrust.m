function [value,isterminal,direction] = event_cont_thrust(t_seconds, x, parameters)

%{ %event function that stops integration when s/c reaches SOI 
%  INPUT : 
% 
% % x : [3,1] distance with respect to the planet
% % data : parameters of the hydraulic plant
% % t : MJD2000
% % % %  OUTPUT 
% pos_diff : distance from SOI
% isterminal : switch to halt integration 
% direction : direction of the states to reach the event
%}

    %parameteres :
%     isEP = parameters.isEP; %1:EP thrust, 0 Chem thust
    t_BO = parameters.t_BO; %burnout time (chemical thrust)
    t0sym = parameters.t0sym*86400; %burnout time (chemical thrust)
    
    isterminal = ones(3,1); %stops integration
   
    muS = astroConstants(4);
    muP = astroConstants(14);
    
    %computation of conditions
    r_sc = kep2car(x(1:6));
    r_sc_norm = norm(r_sc);
    kepPlanet = uplanet(t_seconds/86400,4);
    rP = kep2car2(kepPlanet, muS);
    rP_norm = norm(rP);
    
    
    %SOI definition
    rSOI = rP_norm * (muP/muS)^(2/5);
    

%     %eclipses: 
%     if isEP
%         isterminal(2) = 1;
%         
%         %eclipse check
%         if 
%             eclipse = 1;
%         else
%             eclipse = 0;
%         end
%     else
%         isterminal(2) = 0;
%         eclipse = 1;
%     end
    
    eclipse = 1;
    %check on position 
    value = [rSOI - r_sc_norm; eclipse; (t_seconds - t0sym) - t_BO]; 
    %isterminal defied previously
    direction = [0; 0; 0]; %to be corrected
end

