function [  l_pos_Vect_m, hyp_arrival , kep_hyp_arr] = capture_plot(kep_cap_desired, vv_inf, rM, mu, Thrust , N_firings, parameters)
% function [TT, YY, parout, hyp_arrival , kep_hyp_arr, kep_cap_arr] = capture_plot(kep_cap_desired, vv_inf, rM, mu, Thrust , N_firings, parameters)
    %Plot of asymptotes and hyperbolic path in occurrence of the flyby
    %   INPUT
    %   - delta: deflection angle (free parameters JUST TO RUN SIMULATIONS NOW)
    %   - rM [km] planet position wrt SUN
    %   - kep_cap keplerian parameters of the desired parking orbit
    %OUTPUT
    %   - rp [km] minimum radius of approach of hyperbolic path (from center)
    %   - dt_flyby [hours] TOF across Flyby planet's SOI
    %   - dv_req [km/s] impulsive delta_v
%     delta_t_before_p = parameters.dt_p;
    muS = astroConstants(4);
    %Discretization of hyperbola
    n_iter = 1000;
    %Computation of the hyperbola

    rp = kep_cap_desired(1)*(1 - kep_cap_desired(2));
    e_hyp = 1 + rp*norm(vv_inf)^2/mu;
    h_hyp = sqrt(mu*rp*(1+e_hyp));
    a_hyp =((h_hyp^2)/mu)/(e_hyp^2-1);

    theta_inf = acos(1/e_hyp);
    hyp_arrival.theta_inf = theta_inf;
    hyp_arrival.imp_par = abs(a_hyp) * sqrt(e_hyp^2 - 1);

    %SOI definition
    nrM = norm(rM); 
    rSOI = nrM * (mu/muS)^(2/5);

    correc_hyp_f = @(alpha) correc_SOI_f(rSOI,alpha,a_hyp,e_hyp,theta_inf,mu);
    correc_hyp = fminsearch(correc_hyp_f,0);
    ltheta_hyp = linspace(0,theta_inf+correc_hyp,n_iter);

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

    %Hyperbolas plot
    l_pos_Vect_m = zeros(3,n_iter);

    kep_hyp_arr = kep_cap_desired;
    kep_hyp_arr(1) = -a_hyp;
    kep_hyp_arr(2) = e_hyp;

    [~, v_ell] = kep2car2(kep_cap_desired, mu);
    [~, v_hyp] = kep2car2(kep_hyp_arr, mu);
    dv_req = v_ell - v_hyp;

    %optimal parameters
    hyp_arrival.dv_opt = norm(vv_inf)*sqrt((1-kep_cap_desired(2))/2);
    hyp_arrival.rp_opt = 2 * (1-kep_cap_desired(2))/(1+kep_cap_desired(2)) * mu/norm(vv_inf)^2;
    hyp_arrival.imp_par_opt = rp*sqrt(2/(1 - kep_cap_desired(2)));

    for i = 1:n_iter
        kep_hyp_arr(6) = ltheta_hyp(i);
        [r_m,~]= kep2car2(kep_hyp_arr,mu);
        l_pos_Vect_m(:,i) = r_m;
    end

%     figure
%     hold all
% 
%     r = rSOI;
%     xc = 0;
%     yc = 0;
% 
%     theta = linspace(0,2*pi);
%     x = r*cos(theta) + xc;
%     y = r*sin(theta) + yc;
% 
%     plot(x,y)
%     plot3(l_pos_Vect_m(1,:),l_pos_Vect_m(2,:), l_pos_Vect_m(3,:),'DisplayName','SOI'), hold on
    
%     [X, Y, Z] = ellipsoid(0, 0, 0, r, r, r, 100);    % spheric centered Mars
%     planet = surf(X, Y, -Z,'Edgecolor', 'none','DisplayName','SOI');
%     hold on
%     set(planet,'FaceColor','b')
%     set(planet,'FaceAlpha',0.1)

%     I = imread('Mars.jpg');                            % Mars image
%     RI = imref2d(size(I));
%     RI.XWorldLimits = [-180 180];                       % Mars image x sizes
%     RI.YWorldLimits = [-90 90];                         % Mars image y sizes
%     rM = almanac('Mars','Radius','kilometers','sphere');
%     [X, Y, Z] = ellipsoid(0, 0, 0, rM, rM, rM, 100);    % spheric centered Mars
%     planet = surf(X, Y, -Z,'Edgecolor', 'none','displayname','Mars');
%     hold on
%     set(planet,'FaceColor','texturemap','Cdata',I)
%     axis equal
%     
%     r_desired = zeros(3,360);
%     for jj = 1:360
%         kep_cap_desired(6) = deg2rad(jj);
%         r_desired(:,jj) = kep2car2(kep_cap_desired, mu);
%     end
%     plot3(r_desired(1,:), r_desired(2,:), r_desired(3,:),'g','DisplayName','Desired Arrival Orbit'), hold on
%     plot3(r_desired(1,1), r_desired(2,1), r_desired(3,1), 'o','DisplayName','Hyperbola Pericenter'), hold on
    
    X0_hyp = zeros(7,1);
    [rr0, vv0] = kep2car2(kep_hyp_arr, mu);
    X0_hyp(1:3) = rr0;
    X0_hyp(4:6) = vv0;
    hyp_arrival.kep_hyp_arr = kep_hyp_arr;
    hyp_arrival.dv_req =  norm(dv_req)*1000;  %definition of desired delta_v
    parameters.mP_req = parameters.M0 *(exp(hyp_arrival.dv_req/(9.81*parameters.Isp)) - 1)/ ...
                                        exp(hyp_arrival.dv_req/(9.81*parameters.Isp));
%    
%     parameters.t_BO = parameters.mP_req/norm(Thrust)* parameters.Isp * 9.81;
%     delta_t_before_p = parameters.t_BO/2; 

%     parameters.tmax = parameters.t0sym +  (2*dt_hyp)/86400;
%     X0_hyp(4:6) = -X0_hyp(4:6);
%     X0_hyp = [X0_hyp; 0];      
    %
%     [~, Y] = cart_cont_thrust_model(X0_hyp, parameters);

%     %initial hyperbola plot
%     plot3(Y(:,1), Y(:,2), Y(:,3), ':','DisplayName','Uncontrolled Hyperbolic flight'), hold on

%     parameters.tmax = parameters.t0sym +  (dt_hyp - delta_t_before_p)/86400;                                
%     [T, Y, P] = cart_cont_thrust_model(X0_hyp, parameters);
    
%     TT = T;
%     YY = Y;
%     parout = P;
%     
%     plot3(YY(:,1), YY(:,2), YY(:,3),'DisplayName','Effective Hyperbolic flight'), hold on
%     plot3(YY(1,1), YY(1,2), YY(1,3), 'ro','DisplayName','SOI injection'), hold on

%     % thrust activation
%     parameters.t0sym = parameters.tmax;
%     parameters.tmax = parameters.t0sym + parameters.t_BO/86400;
%     parameters.T = Thrust;              % thrust [N] (constant for the moment)
% 
%     T0_firing = parameters.t0sym*86400;
%     Tend_firing = parameters.tmax*86400;
% 
%     [T, Y, P] = cart_cont_thrust_model(YY(end,:), parameters);
% 
%     plot3(Y(:,1), Y(:,2), Y(:,3), 'r','DisplayName',strcat('Firing (T = ', num2str(norm(parameters.T)), 'N, tBO = ', num2str(norm(parameters.t_BO)),'s)')), hold on
%     plot3(Y(1,1), Y(1,2), Y(1,3), 'o', 'MarkerSize',15,'DisplayName','Firing start'), hold on
%     plot3(Y(end,1), Y(end,2), Y(end,3), 'ro', 'MarkerSize',15,'DisplayName','Firing end'), hold on

%     YY = [YY; Y];
%     TT = [TT; T];
%     parout = [parout; P];

    % Thrust de-activated

%     kep_cap_arr = car2kep(YY(end,1:3), YY(end,4:6), mu);
%     T_orb = 2*pi*sqrt(kep_cap_arr(1)^3/mu);

%     parameters.t0sym = parameters.tmax;
%     parameters.tmax = parameters.t0sym + T_orb/86400;
%     parameters.T = [0; 0; 0];            

%     [T, Y_arr, P] = cart_cont_thrust_model(YY(end,:), parameters);

%     YY = [YY; Y_arr];
%     TT = [TT; T];
%     parout = [parout; P];

%     figure1 = figure(1);
%     plot3(Y_arr(:,1),Y_arr(:,2),Y_arr(:,3),'b','DisplayName','Effective Capture Ellipse'), hold on
%     legend()
% 
%     annotation(figure1,'textbox',...
%     [0.075 0.65 0.15 0.3],...
%     'String', ...
%     {'Desired Orbit Kep:', ...
%     strcat('a: ', num2str(kep_cap_desired(1)), 'km'), ...
%     strcat('e: ', num2str(kep_cap_desired(2))), ...
%     strcat('i: ', num2str(rad2deg(kep_cap_desired(3))), '�'), ...
%     strcat('OM: ', num2str(rad2deg(kep_cap_desired(4))), '�'), ...
%     strcat('om: ', num2str(rad2deg(kep_cap_desired(5))), '�'), ...
%     strcat('theta: ', num2str(rad2deg(kep_cap_desired(6))), '�'), ...
%     ' ', ...
%     'Effective Orbit Kep:', ...
%     strcat('a: ', num2str(kep_cap_arr(1)), 'km'), ...
%     strcat('e: ', num2str(kep_cap_arr(2))), ...
%     strcat('i: ', num2str(rad2deg(kep_cap_arr(3))), '�'), ...
%     strcat('OM: ', num2str(rad2deg(kep_cap_arr(4))), '�'), ...
%     strcat('om: ', num2str(rad2deg(kep_cap_arr(5))), '�'), ...
%     strcat('theta: ', num2str(rad2deg(kep_cap_arr(6))), '�'), ...
%     ' ', ...
%     'Engine Specifics:', ...
%      strcat('Isp',num2str(parameters.Isp), 's'), ...
%      strcat('T',num2str(Thrust),' N'), ...
%     'N:firings:',num2str(N_firings), ...
%     ' ' , ...
%     'SM:', ...
%     strcat('M0', num2str(parameters.M0), 'kg'), ...
%     strcat('Mp damped', num2str(YY(end,7)), 'kg')} , ...
%     'FitBoxToText','on');
% 
% 
%     %delta v
%     figure()
%     delta_v_eff = parameters.Isp * 9.81 * log(parameters.M0./(parameters.M0 - YY(:,7)));
%     plot(TT, delta_v_eff), hold on, title('Delta v profile'), ylabel('Delta v [m/s]'), xlabel('t [s]')
%     yline(hyp_arrival.dv_req ,'b-'), hold on,
%     xline(T0_firing,'r'), hold on, xline(Tend_firing,'r'), hold off
%     legend('Delta_v(t)', 'Required Impulsive Delta v', strcat('t_BO:', num2str(parameters.t_BO) ,' s'))
%     %keplerian parameters
%     figure()
%     sgtitle('Keplerian Parameters variation')
%     ttx{1} = 'a'; ttx{2} = 'e';
%     ttx{3} = 'i'; ttx{4} = 'OM'; 
%     ttx{5} = 'om'; ttx{6} = 'th';

%     YY_kep = YY;
%     for kk = 1:length(TT)
%         YY_kep(kk,1:6) = car2kep(YY(kk,1:3), YY(kk,4:6), mu);
%     end
%     YY_kep(:,6) = real(YY_kep(:,6));
%     for kk = 1:6
%         YY_kep(end,kk)
%         subplot(2,3,kk), plot(TT, YY_kep(:,kk));
%         title (ttx{kk})
%     end
end

% Function used to get a correction on the theta for the hyperbola arc to
% reach the Sphere of Influence
function val = correc_SOI_f(rSOI,alpha,a,e,theta,muM)
    [r,~]= kep2car2([a,e,0,0,0,theta+alpha],muM);
    val = (norm(r)-rSOI)^2;
end



