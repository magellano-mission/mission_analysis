function [T_tot,states_tot, parout_tot, events] = model_continuous_thrust_gauss(X0, parameters)
        isTactive = 1;
        T_tot = []; states_tot = []; parout_tot = []; events = [];
        T = parameters.t0sym*86400;
        
        while T(end) <  parameters.tmax*86400 || i_e == 1
                        
            [T, states, parout, i_e] = gauss_cont_thrust_model(X0, parameters);
            
            X0 = states(end,:); 
            parameters.t0sym = T(end)/86400; 
            parameters.M0 = parameters.M0 - isTactive*states(end,7);
            
            if i_e == 3 %t_BO is over (add mP is over)
                parameters.T = [0; 0; 0];
                parameters.t_BO = parameters.tmax*86400 - T(end) + 1;
                isTactive = 0;
            else
                parameters.t_BO = parameters.t_BO - (T(end) - parameters.t0sym);
            end
%             if i_e == 2 %eclipse presence -> EP turned off (it must be
%             turned onanyway, will code it later)
%                 parameters.T = [0; 0; 0];
%             end
            if i_e == 4
                parameters.T = [0; 0; 0];
                parameters.t_BO = parameters.tmax*86400 - T(end) + 1;
                isTactive = 0;
            end
            T_tot = [T_tot; T];
            states_tot = [states_tot; states];
            parout_tot = [parout_tot; parout];
            events = [T(end) i_e];
        end


end

