function [T_tot, states_tot, parout_tot, events] = model_continuous_thrust_cart(X0, parameters)
        %{
        %simulation starts with Thrust active
        %  legend of events (i_e)
        % 1: SOI reached;
        % 2: under eclipse cone
        % 3: Burn Out time reached
        % 4: Delta_v required is over
        %}
        %initialization of all needed vectors
        isTactive = 1; i_e = 0;
        T_tot = []; states_tot = []; parout_tot = []; events = [];
        T = parameters.t0sym*86400;
        
        while T(end) <  parameters.tmax*86400 || i_e == 1
                        
            [T, states, parout, i_e] = cart_cont_thrust_model(X0, parameters);
            
            %increment all the vectors
            T_tot = [T_tot; T];
            states_tot = [states_tot; states];
            parout_tot = [parout_tot; parout];
            events = [events; T(end) i_e]; %keeps track of events and time of events
            
            %inizialization of next iteration
            X0 = states(end,:); 
            parameters.t0sym = T(end)/86400;  
            
            if isTactive
               i_e = 0;
%                ind_eclipse = 0; %no eclipse coming
%           if i_e == 2 %eclipse presence -> EP turned off (it must be
%             turned onanyway, will code it later)
%             parameters.T = [0; 0; 0];
%             ind_eclipse = 1;  
%           end
            
            if i_e == 3                %  t_BO is over (therefore mP is over)
                parameters.T = [0; 0; 0];
                parameters.t_BO = Inf; %parameters.tmax*86400 - T(end) + 1;%prevents from stopping again because of t_BO
                isTactive = 0;
            else
                %reduction of t_BO
                parameters.t_BO = parameters.t_BO - (T(end) - parameters.t0sym);
            end

            if i_e == 4 %delta_v required is reached
                parameters.T = [0; 0; 0];
                parameters.t_BO = parameters.tmax*86400 - T(end) + 1;
                isTactive = 0;
            end
          %   when there are eclipses and you come out of cone of eclipse
%             else %(isTactive false)
%                 if ind_eclipse == 1 && i_e ~= 2
%                  isTactive = 1;
%                  parameters.T = [0; 0; 0];
%                 end
%             end
        end


end

