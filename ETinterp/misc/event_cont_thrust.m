function [value,isterminal,direction] = event_cont_thrust(t_seconds, x, parameters)

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
%     isEP = parameters.isEP; %1:EP thrust, 0 Chem thust
    t_BO = parameters.t_BO; %burnout time (chemical thrust)
    t0sym = parameters.t0sym*86400; %burnout time (chemical thrust)
    
    isterminal = ones(4,1); %stops integration if 1
   
    muS = astroConstants(4);
    muP = astroConstants(14);
    
    %computation of conditions
    r_sc = x(1:3);              %CARTESIAN! 
    r_sc_norm = norm(r_sc);
    kepPlanet = uplanet(t_seconds/86400,4);
    rP = kep2car2(kepPlanet, muS);
    rP_norm = norm(rP);

    %SOI definition
    rSOI = rP_norm * (muP/muS)^(2/5);
       
    delta_v_actual = parameters.Isp * 9.81* log(parameters.M0/(parameters.M0 - x(7))); %[m/s]
    %checks
    value = [rSOI - r_sc_norm; ...
             (t_seconds - t0sym) - t_BO; ...
             delta_v_actual - parameters.delta_v_req]; 
         
    %isterminal defied previously
    direction = [0; 0; 0];
end

