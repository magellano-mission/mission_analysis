function [T, Y, parout, t_e, y_e, i_e] = cart_cont_thrust_model(X0, parameters)
 function [dxdt, parout] = cart_cont_thrust(t_seconds, x, parameters)

    %parameteres :
    Thrust = parameters.T;           %thrust [N] (constant for the moment)
    Isp = parameters.Isp;              %specific impulse [s]
    M0 = parameters.M0;              %structural mass
    isInterp = parameters.isInterp;
    
    dxdt = zeros(length(x),1);

    if isInterp == 1
        mu = astroConstants(4);
    else
        mu = astroConstants(14);
        J2 = 0.00196; %Mars J2
    end
    
    in_eclipse = 0;
    
    %States definition
    rr = x(1:3); if size(rr,1) < size(rr,2) rr = rr'; end
    vv = x(4:6); if size(vv,1) < size(vv,2) vv = vv'; end
    %Damped prop mass (>0)
    Mprop = x(7);
    
    % Rotation Matrix 
    t_hat = vv/norm(vv);
    h_hat = cross(rr,vv)/norm(cross(rr,vv));
    n_hat = cross(h_hat,t_hat);
    A = [ t_hat  n_hat  h_hat ];
    
    % Thrust Vector (add burnout time - callback and eclipses)
    M = M0 - Mprop;
    a_thrust = ((in_eclipse == 0) && (parameters.isEP == 0)) * Thrust/M/1000;

    %DISTURBANCES
    if parameters.isPerturbed
    if ~isInterp
        t0sym = parameters.t0sym;
        % Compute altitude from ellipsoid model 
        [~, R_real] = MCI2altitude(rr);    %[km]
        r = norm(x(1:3));
        % J2 Effect 
        k = 3/2 .* J2 .* mu .* R_real^2/r^4;
        a_j2_x = k*(rr(1)/r .*(5 *rr(3)^2/r^2 -1));
        a_j2_y = k*(rr(2)/r .*(5 *rr(3)^2/r^2 -1));
        a_j2_z = k*(rr(3)/r .*(5 *rr(3)^2/r^2 -3));
        a_J2 = [a_j2_x; a_j2_y; a_j2_z]; % [km/s^2]

        %OTHER GRAVITY EFFECTS
        %Conversion day to seconds
        %day2sec = 86400;

        % PHOBOS
        M_phobos = 1.07e16;
        G = astroConstants(1);
        mu_phobos = M_phobos .* G;

        % Retrieving moon position
        %mjd = date2mjd2000(data.date) + t/day2sec;
        stat_phobos = [9376, 0.0151, deg2rad(1.075) deg2rad(317.671) deg2rad(150.057) 0];
        [~, x_phobos] = ode113(@gauss_sat, [t0sym*86400 t_seconds+1],  stat_phobos);
        [rr_ph, ~] = kep2car2(x_phobos(end,:),mu);
        r_ph = norm(rr_ph); %[km]

        % Retrieving relative position moon / spacecraft
        rr_rel_ph = rr_ph - rr;
        r_rel_ph = norm(rr_rel_ph);

        % Moon gravitational acceleration
        a_phobos = mu_phobos .* (rr_rel_ph/r_rel_ph^3 - rr_ph/r_ph^3);

        % DEIMOS
        M_deimos = 1.4762e15;
        G = astroConstants(1);
        mu_deimos = M_deimos .* G;

        % Retrieving moon position
        %mjd = date2mjd2000(data.date) + t/day2sec;
        stat_deimos = [23458, 0.0002, deg2rad(1.788) deg2rad(316.657) deg2rad(260.729) 0];
        [~, x_deimos] = ode113(@gauss_sat, [t0sym*86400 t_seconds+1],  stat_deimos);
        [rr_d, ~]=kep2car2(x_deimos(end,:), mu);
        r_d = norm(rr_d); %[km]

        % Retrieving relative position moon / spacecraft
        rr_rel_d = rr_d - rr;
        r_rel_d = norm(rr_rel_d);

        % Moon gravitational acceleration
        a_deimos = mu_deimos .* (rr_rel_d/r_rel_d^3 - rr_d/r_d^3);
        % %ADD ATMOSPHERE MODEL (might be unuseful)
        % ADD SRP 
        
         %Total of accelerations
        a_tot = a_J2 + a_phobos + a_deimos;
    else
        % JUPITER
        mu_J = astroConstants(15);
        kep_J = uplanet(t_seconds/86400,5);
        [rr_J, ~] = kep2car2(kep_J, mu);
        r_J = norm(rr_J); %[km]

        % Retrieving relative position moon / spacecraft
        rr_rel_J = rr_J - rr;
        r_rel_J = norm(rr_rel_J);

        % Moon gravitational acceleration
        a_Jupiter = mu_J .* (rr_rel_J/r_rel_J^3 - rr_J/r_J^3);   
        
        % EARTH
        mu_E = astroConstants(13);
        kep_E = uplanet(t_seconds/86400, 5);
        [rr_E, ~] = kep2car2(kep_E, mu);
        r_E = norm(rr_E); %[km]

        % Retrieving relative position moon / spacecraft
        rr_rel_E = rr_E - rr;
        r_rel_E = norm(rr_rel_E);

        % Moon gravitational acceleration
        a_Earth = mu_E .* (rr_rel_E/r_rel_E^3 - rr_E/r_E^3);   
        
        %SRP 
        a_SRP = solarPressure(rr, M, parameters);
        
        a_tot = a_Jupiter + a_Earth + a_SRP;
    end
    else
        a_tot = zeros(3,1);
    end
    a_thrust_cart = A*a_thrust;
    a_tot = a_tot + a_thrust_cart;
    % Derivative of the keplerian parameters  
    dxdt(1) = x(4);
    dxdt(2) = x(5);
    dxdt(3) = x(6);
    dxdt(4) = -mu/sqrt(sum(x(1:3).^2))^3*x(1) + a_tot(1);
    dxdt(5) = -mu/sqrt(sum(x(1:3).^2))^3*x(2) + a_tot(2);
    dxdt(6) = -mu/sqrt(sum(x(1:3).^2))^3*x(3) + a_tot(3);
    dxdt(7) =  norm(Thrust)/Isp/9.81;
    
    delta_v = Isp* 9.81 .* log(M0/ (M0 - x(7)));
%     Itot = Itot + trapz([tlast t_seconds],[norm(Thrust) norm(Thrust)]);
    parout = [M Itot delta_v];

 end


%%%% Callback nested 
    function state = callback( t , x, flag, ~ )
        if isempty(flag)
            tlast = t(end);
        end
        state = 0 ;
    end

tlast = 0;                              % integration time
Itot = 0;                           % integral over the time of Vout

if parameters.event == 1
    [T, Y, t_e, y_e, i_e] = ode15s(@cart_cont_thrust,[parameters.t0sym parameters.tmax]*86400, X0, parameters.opt, parameters);
else
    [T, Y] = ode15s(@cart_cont_thrust,[parameters.t0sym parameters.tmax]*86400, X0, parameters.opt, parameters);
end
parout = [];
for kkk=1:length(T)
    [~, P] = cart_cont_thrust(T(kkk), Y(kkk,:), parameters);
    parout = [parout; P];
end

end
