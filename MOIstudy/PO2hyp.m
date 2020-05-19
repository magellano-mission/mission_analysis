function [ output_pos, hyp_arrival , kep_hyp_dep] = PO2hyp(kep_departure, vv_inf, rM, mu, parameters, Thrust, N_firings, flag_maneuver, plot_flag)
% function [TT, YY, parout, hyp_arrival , kep_hyp_arr, kep_cap_arr] = capture_plot(kep_cap_desired, vv_inf, rM, mu, Thrust , N_firings, parameters)
    %Plot of asymptotes and hyperbolic path in occurrence of the flyby
    %   INPUT
    %   - kep_cap_desired : keplerian states vector for the orbit aimed
    %   - vv_inf: [3x1] vector of planetocentric approach velocity
    %   - rM [km] helyocentric planet position
    %   - mu  planet gravity constant
    %   - parameters....
    %   - Thrust, N_firings ...
    %   - flag_maneuver: 'arrival' or [] (departure)
    %   - plot_flag: 0:no plot, 1:geometrical plot, 2:integration
    %   OUTPUT
    %   - output_pos: states (plot_flag == 0->analytical, plot_flag == 1-> numerical)
    %   - hyp_arrival
    %   - dv_req [km/s] impulsive delta_v

    if nargin == 8
        plot_flag = 0;
    end

    muS = astroConstants(4);
    %Discretization of hyperbola
    n_iter = 1000;
    %Computation of the hyperbola
    
    rp = kep_departure(1)*(1 - kep_departure(2));
    e_hyp = 1 + rp*norm(vv_inf)^2/mu;
    h_hyp = sqrt(mu*rp*(1+e_hyp));
    a_hyp =((h_hyp^2)/mu)/(e_hyp^2-1);
    
    theta_inf = acos(1/e_hyp);
    hyp_arrival.theta_inf = theta_inf;
    hyp_arrival.imp_par = abs(a_hyp) * sqrt(e_hyp^2 - 1);
  
    %SOI definition
    nrM = norm(rM); 
    rSOI = nrM * (mu/muS)^(2/5);

    %Hyperbolas plot
    l_pos_Vect_m = zeros(3,n_iter);
    
    kep_departure(6) = 0;
    kep_hyp_dep = kep_departure;
    kep_hyp_dep(1) = a_hyp; %>0
    kep_hyp_dep(2) = e_hyp; %>1
    
    correc_hyp = acos((-a_hyp*(1-e_hyp^2)/rSOI - 1)/e_hyp);

%     correc_hyp_f = @(alpha) correc_SOI_f(rSOI,alpha, kep_hyp_arr, theta_inf, mu);
%     correc_hyp = fminsearch(correc_hyp_f,0); 


    I = imread('Earth.jpg');                            
    rP = almanac('Earth','Radius','kilometers','sphere');
    ltheta_hyp = linspace(0,correc_hyp,n_iter);
    
    lt_hyp = zeros(1,n_iter);
    
    %TOF computation
    for i = 1:n_iter
        %eccentric anomalies
        f_hyp = 2*atanh(sqrt((e_hyp-1)/(e_hyp+1))*tan(ltheta_hyp(i)/2));
        m_hyp = e_hyp*sinh(f_hyp) - f_hyp;
        lt_hyp(i) = m_hyp*(h_hyp^3/mu^2)/((e_hyp^2-1)^(3/2));
    end
    
    %Computation of TOF
    dt_hyp = abs(lt_hyp(1) - lt_hyp(n_iter)); % [s]
    hyp_arrival.dt_hyp = dt_hyp;
    
    [r0_dep, v_ell] = kep2car2(kep_departure, mu);
    [~, v_hyp] = kep2car2(kep_hyp_dep, mu);
    dv_req = v_hyp - v_ell;

    %optimal parameters (only for capture!)
    hyp_arrival.dv_opt = norm(vv_inf)*sqrt((1-kep_departure(2))/2);
    hyp_arrival.rp_opt = 2 * (1-kep_departure(2))/(1+kep_departure(2)) * mu/norm(vv_inf)^2;
    hyp_arrival.imp_par_opt = rp*sqrt(2/(1 - kep_departure(2)));

    for i = 1:n_iter
        kep_hyp_dep(6) = ltheta_hyp(i);
        [r_m,~]= kep2car2(kep_hyp_dep,mu);
        l_pos_Vect_m(:,i) = r_m;
    end
    output_pos = l_pos_Vect_m';
    
    kep_hyp_dep(6) = ltheta_hyp(1);

    hyp_arrival.kep_hyp_arr = kep_hyp_dep;
    hyp_arrival.dv_req =  norm(dv_req)*1000;  %definition of desired delta_v
    parameters.mP_req = parameters.M0 *(exp(hyp_arrival.dv_req/(9.81*parameters.Isp)) - 1)/ ...
                                        exp(hyp_arrival.dv_req/(9.81*parameters.Isp));

%   preliminary computation of required propellant
    parameters.t_BO = parameters.mP_req/norm(Thrust)* parameters.Isp * 9.81; %computation of rectangular profile Thrust 
if plot_flag %plot of geometrical orbits
    figure
    xc = 0;
    yc = 0;

    theta = linspace(0,2*pi);
    x = rSOI*cos(theta) + xc;
    y = rSOI*sin(theta) + yc;

%   equatorial plane
    plot(x,y, 'DisplayName','SOI boundary'), 
    plot3(l_pos_Vect_m(1,:),l_pos_Vect_m(2,:), l_pos_Vect_m(3,:),'DisplayName','Analytical  Hyperbola'), hold on
%     
    [X, Y, Z] = ellipsoid(0, 0, 0, rSOI, rSOI, rSOI, 100);    % spheric centered Mars
    planet = surf(X, Y, -Z,'Edgecolor', 'none','DisplayName','SOI');
    hold on
    set(planet,'FaceColor','b')
    set(planet,'FaceAlpha',0.1)
    
%   Planet plot
    RI = imref2d(size(I));
    RI.XWorldLimits = [-180 180];                       % Mars image x sizes
    RI.YWorldLimits = [-90 90];                         % Mars image y sizes
    [X, Y, Z] = ellipsoid(0, 0, 0, rP, rP, rP, 100);    % spheric centered Mars
    planet = surf(X, Y, -Z,'Edgecolor', 'none','displayname','Mars');
    hold on
    set(planet,'FaceColor','texturemap','Cdata',I)
    axis equal
    
%   desired capture orbit plot
    r_departure = zeros(3,360);
    vvvv = r_departure;
    for jj = 1:360
        kep_departure(6) = deg2rad(jj-1);
        [r_departure(:,jj), vvvv(:,jj)] = kep2car2(kep_departure, mu);
    end
    kep_departure(6) = 0;
    plot3(r_departure(1,:), r_departure(2,:), r_departure(3,:),'g','DisplayName','Desired Arrival Orbit'), hold on
    plot3(r_departure(1,1), r_departure(2,1), r_departure(3,1), 'o','DisplayName','PO Pericenter'), hold on
    xlim([-13000 13000]), ylim([-13000 13000]), zlim([-13000 13000])
    if plot_flag == 2
%   integration of hyperbolic fligt
    X0_PO = zeros(7,1);
    X0_PO(1:3) = r0_dep;
    X0_PO(4:6) = v_ell;

    delta_t_before_p = parameters.perc_tBO_before_p * parameters.t_BO; %how much time before reaching the pericenter occurs the burn
    
    T_orb = 2*pi*sqrt(kep_departure(1)^3/mu);
%   integration of complete hyperbola
    parameters.tmax = parameters.t0sym +  (T_orb)/86400; %[MJD2000]
%     X0_hyp(4:6) = -X0_hyp(4:6); 
    parameters.T = zeros(3,1);  
    [~, tot_departure] = cart_cont_thrust_model(X0_PO, parameters);

%   complete PO plot
    plot3(tot_departure(:,1), tot_departure(:,2), tot_departure(:,3), 'k','DisplayName','Uncontrolled flight @PO'), hold on

%   no-thrust path integration
    parameters.tmax = parameters.t0sym +  (T_orb - delta_t_before_p)/86400;                                
    [T_no_control, Y_no_control, P_no_control] = cart_cont_thrust_model(X0_PO, parameters);
    
    TT = T_no_control;
    YY = Y_no_control;
    parout = P_no_control;
     
    plot3(YY(:,1), YY(:,2), YY(:,3), 'DisplayName','Coasting maneuver'), hold on
    plot3(YY(1,1), YY(1,2), YY(1,3), 'r*','DisplayName','Departure'), hold on

%     thrust activation
    parameters.t0sym = parameters.tmax;                             %[MJD200]
    parameters.tmax = parameters.t0sym + parameters.t_BO/86400;     %[MJD200]
    parameters.T = Thrust;              % thrust [N] (constant profile for the moment)

    T0_firing = parameters.t0sym*86400;                             %[s]
    Tend_firing = parameters.tmax*86400;                            %[s]

    [T_thrust, Y_thrust, P_thrust] = cart_cont_thrust_model(YY(end,:), parameters);

    plot3(Y_thrust(:,1), Y_thrust(:,2), Y_thrust(:,3), 'r','DisplayName',strcat('Firing (T = ', num2str(norm(parameters.T)), 'N, tBO = ', num2str(norm(parameters.t_BO)),'s)')), hold on
    plot3(Y_thrust(1,1), Y_thrust(1,2), Y_thrust(1,3), 'o', 'MarkerSize',15,'DisplayName','Firing start'), hold on
%     plot3(Y_thrust(end,1), Y_thrust(end,2), Y_thrust(end,3), 'ro', 'MarkerSize',15,'DisplayName','Firing end'), hold on

    YY = [YY; Y_thrust];
    TT = [TT; T_thrust];
    parout = [parout; P_thrust];

%   Hyperbolic Path
    kep_cap_arr = car2kep(YY(end,1:3), YY(end,4:6), mu);

    parameters.t0sym = parameters.tmax;
    parameters.tmax = parameters.t0sym + dt_hyp/86400;
    parameters.T = [0; 0; 0];  

    [T_hyp, Y_hyp, P_hyp] = cart_cont_thrust_model(YY(end,:), parameters);

    YY = [YY; Y_hyp];
    TT = [TT; T_hyp];
    parout = [parout; P_hyp];

%   plot of capture orbit result
    figure1 = gcf;
    plot3(Y_hyp(:,1),Y_hyp(:,2),Y_hyp(:,3),'b','DisplayName','Effective Departure Hyperbola'), hold on
%     xlim([-11000 11000]), ylim([-11000 11000]), zlim([-11000 11000])
    xlim([-rSOI rSOI]), ylim([-rSOI rSOI]), zlim([-rSOI rSOI])

    legend()

    annotation(figure1,'textbox',...
    [0.075 0.65 0.15 0.3],...
    'String', ...
    {'Desired Orbit Kep:', ...
    strcat('a: ', num2str(kep_departure(1)), 'km'), ...
    strcat('e: ', num2str(kep_departure(2))), ...
    strcat('i: ', num2str(rad2deg(kep_departure(3))), '°'), ...
    strcat('OM: ', num2str(rad2deg(kep_departure(4))), '°'), ...
    strcat('om: ', num2str(rad2deg(kep_departure(5))), '°'), ...
    strcat('theta: ', num2str(rad2deg(kep_departure(6))), '°'), ...
    ' ', ...
    'Effective Orbit Kep:', ...
    strcat('a: ', num2str(kep_cap_arr(1)), 'km'), ...
    strcat('e: ', num2str(kep_cap_arr(2))), ...
    strcat('i: ', num2str(rad2deg(kep_cap_arr(3))), '°'), ...
    strcat('OM: ', num2str(rad2deg(kep_cap_arr(4))), '°'), ...
    strcat('om: ', num2str(rad2deg(kep_cap_arr(5))), '°'), ...
    strcat('theta: ', num2str(rad2deg(kep_cap_arr(6))), '°'), ...
    ' ', ...
    'Engine Specifics:', ...
     strcat('Isp',num2str(parameters.Isp), 's'), ...
     strcat('T',num2str(Thrust),' N'), ...
    'N:firings:',num2str(N_firings), ...
    ' ' , ...
    'SM:', ...
    strcat('M0', num2str(parameters.M0), 'kg'), ...
    strcat('Mp damped', num2str(YY(end,7)), 'kg')} , ...
    'FitBoxToText','on');


%   delta v
    figure()
    delta_v_eff = parameters.Isp * 9.81 * log(parameters.M0./(parameters.M0 - YY(:,7)));
    plot(TT, delta_v_eff, 's'), hold on, title('Delta v history'), ylabel('Delta v [m/s]'), xlabel('t [s]')
    yline(hyp_arrival.dv_req ,'b-'), hold on,
    xline(T0_firing,'r'), hold on, xline(Tend_firing,'r'), hold off
    legend('Delta_v(t)', 'Required Impulsive Delta v', strcat('t_BO:', num2str(parameters.t_BO) ,' s'))

%   keplerian parameters history
    figure()
    sgtitle('Keplerian Parameters variation')
    ttx{1} = 'a'; ttx{2} = 'e';
    ttx{3} = 'i'; ttx{4} = 'OM'; 
    ttx{5} = 'om'; ttx{6} = 'th';

    YY_kep = YY;
    for kk = 1:length(TT)
        YY_kep(kk,1:6) = car2kep(YY(kk,1:3), YY(kk,4:6), mu);
    end
    YY_kep(:,6) = real(YY_kep(:,6));
    for kk = 1:6
        YY_kep(end,kk)
        subplot(2,3,kk), plot(TT, YY_kep(:,kk));
        title (ttx{kk})
    end
    output_pos = YY;
    end
end
end



