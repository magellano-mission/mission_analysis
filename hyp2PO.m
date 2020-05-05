function [TT, YY, dt_hyp, dv_req, theta_inf, kep_hyp_arr, kep_cap_arr] = hyp2PO(kep_cap_desired, delta,rM, mu, Thrust , parameters)
    %Hyperbolic path in occurrence of the flyby
    %   INPUT
    %   - delta: deflection angle (free parameter JUST TO RUN SIMULATIONS NOW)
    %   - rM [km] planet position wrt SUN
    %   - kep_cap keplerian parameters of the desired parking orbit
    %OUTPUT
    %   - rp [km] minimum radius of approach of hyperbolic path (from center)
    %   - dt_flyby [hours] TOF across Flyby planet's SOI
    %   - dv_req [km/s] impulsive delta_v
    delta_t_before_p = parameters.dt_p;
    muS = astroConstants(4);
    %Discretization of hyperbola
    n_iter = 1000;

    %Computation of the hyperbola
    rp = kep_cap_desired(1)*(1 - kep_cap_desired(2));
    e_m = 1/sin(delta/2);
    h_m = sqrt(mu*rp*(1+e_m));
    a_m =((h_m^2)/mu)/(e_m^2-1);
    theta_inf = acos(1/e_m);

    nrM = norm(rM); 
    %SOI definition
    rSOI = nrM * (mu/muS)^(2/5);
    correc_m_f = @(alpha) correc_SOI_f(rSOI,alpha,a_m,e_m,theta_inf,mu);
    correc_m = fminsearch(correc_m_f,0);
    ltheta_m = linspace(0,theta_inf+correc_m,n_iter);
    lt_m = zeros(1,n_iter);
    
    %TOF computation
    for i = 1:n_iter
        %eccentric anomalies
        f_m = 2*atanh(sqrt((e_m-1)/(e_m+1))*tan(ltheta_m(i)/2));
        m_m = e_m*sinh(f_m) - f_m;
        lt_m(i) = m_m*(h_m^3/mu^2)/((e_m^2-1)^(3/2));
    end

    %Computation of TOF  %[s]
    dt_hyp = abs(lt_m(1) - lt_m(n_iter)); 

    %Hyperbolas plot
    l_pos_Vect_m = zeros(3,n_iter);

    kep_hyp_arr = kep_cap_desired;
    kep_hyp_arr(1) = -a_m;
    kep_hyp_arr(2) = e_m;

    [~, v_ell] = kep2car2(kep_cap_desired, mu);
    [~, v_hyp] = kep2car2(kep_hyp_arr, mu);
    dv_req = v_ell - v_hyp;

    for i = 1:n_iter
        kep_hyp_arr(6) = ltheta_m(i);
        [r_m,~]= kep2car2(kep_hyp_arr,mu);
        l_pos_Vect_m(:,i) = r_m;
    end
         
    X0_hyp = zeros(7,1);
    [rr0, vv0] = kep2car2(kep_hyp_arr, mu);
    X0_hyp(1:3) = rr0;
    X0_hyp(4:6) = vv0;
    
    %definition of desired delta_v
    parameters.delta_v_req =  norm(dv_req)*1000;       % [m/s]
    X0_hyp(4:6) = -X0_hyp(4:6); %%%% it has to disappear
    X0_hyp = [X0_hyp; 0];      
%
    parameters.tmax = parameters.t0sym +  (dt_hyp - delta_t_before_p)/86400;                                
    [T, Y] = cart_cont_thrust_model(X0_hyp, parameters);
    TT = T;
    YY = Y;

    %%%%%%%%%%%%% thrust activation (possible looping with function under developement)
    parameters.t0sym = parameters.tmax;
    parameters.tmax = parameters.t0sym + parameters.t_BO/86400;
    parameters.T = Thrust;              % thrust [N] (constant for the moment)

    [T, Y] = cart_cont_thrust_model(YY(end,:), parameters);

    YY = [YY; Y];
    TT = [TT; T];

    %%%%%%%%%%%% Thrust de-activated

    kep_cap_arr = car2kep(YY(end,1:3), YY(end,4:6), mu);
    T_orb = 2*pi*sqrt(kep_cap_arr(1)^3/mu);
    
    parameters.t0sym = parameters.tmax;
    parameters.tmax = parameters.t0sym + T_orb/86400;
    parameters.T = [0; 0; 0];            

    [T, Y_arr] = cart_cont_thrust_model(YY(end,:), parameters);

    YY = [YY; Y_arr];
    TT = [TT; T];

end

% Function used to get a correction on the theta for the hyperbola arc to
% reach the Sphere of Influence
function val = correc_SOI_f(rSOI,alpha,a,e,theta,muM)
    [r,~]= kep2car2([a,e,0,0,0,theta+alpha],muM);
    val = (norm(r)-rSOI)^2;
end



