%% CAPTURE STUDY
%Propagation of a single-burn maneuver for a generic capture (more refined
%model coming in one day)
%change the Thrust vector, M0, t_BO, Isp, (Mp0 is an output at the moment,
%myabe it is better to understand the order of magnitudes involved)
%delta_t_before_pericenter: trying to understand the impact of the firing
%timing

close all, clear, clc

parameters.isInterp = 0;
parameters.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
parameters.isPerturbed = 0;                                % perturbation switch 0 ,too slow :(
parameters.t_BO = 30*60;                                   % burnout time (chemical thrust)
parameters.T = [0; 0; 0];                                  % thrust [N] (@TNH) (constant profile for the moment)
Thrust = [-3200; 0; 0];                                    % thrust [N] (@TNH) (constant profile for the moment)
parameters.Isp = 280;                                      % specific impulse [s]
parameters.M0 = 5000;                                      % Total Mass of the s/c [kg]
parameters.c_r = 0.5;
parameters.Across_sun = 10;                                % Cross area related to the sun [m^2]
parameters.t0sym = date2mjd2000([2021, 1, 1, 0, 0, 0]);    
delta_t_before_p = parameters.t_BO/2; 
parameters.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12);   %,'Events', @event_cont_thrust);

%capture delta_v
kep_cap = [15000 0.1 deg2rad(15) 0 0 0];
delta = deg2rad(45);
[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);
mu = astroConstants(14);

[dt_hyp, dv_req, theta_inf, kep_hyp] = capture(kep_cap, delta,rM, mu);
%the idea is to start the burn at (t_pericenter - t_BO/2) to emulate at
%most an impulsive maneuver at pericenter

X0_hyp = zeros(7,1);
[rr0, vv0] = kep2car2(kep_hyp, mu);
X0_hyp(1:3) = rr0;
X0_hyp(4:6) = vv0;

parameters.delta_v_req =  norm(dv_req);  %definition of desired delta_v
parameters.tmax = parameters.t0sym +  (2*dt_hyp)/86400;
X0_hyp(4:6) = -X0_hyp(4:6);
X0_hyp = [X0_hyp; 0];      
%
[~, Y] = cart_cont_thrust_model(X0_hyp, parameters);

%initial hyperbola plot
plot3(Y(:,1), Y(:,2), Y(:,3), ':','DisplayName','Uncontrolled Hyperbolic flight'), hold on

parameters.tmax = parameters.t0sym +  (dt_hyp - delta_t_before_p)/86400;                                
[T, Y] = cart_cont_thrust_model(X0_hyp, parameters);
TT = T;
YY = Y;

plot3(YY(:,1), YY(:,2), YY(:,3),'DisplayName','Effective Hyperbolic flight'), hold on
plot3(YY(1,1), YY(1,2), YY(1,3), 'ro','DisplayName','SOI injection'), hold on

%% thrust activation
parameters.t0sym = parameters.tmax;
parameters.tmax = parameters.t0sym + parameters.t_BO/86400;
parameters.T = Thrust;              % thrust [N] (constant for the moment)

[T, Y] = cart_cont_thrust_model(YY(end,:), parameters);

plot3(Y(:,1), Y(:,2), Y(:,3), 'r','DisplayName',strcat('Firing (T = ', num2str(norm(parameters.T)), 'N, tBO = ', num2str(norm(parameters.t_BO)),'s)')), hold on
plot3(Y(1,1), Y(1,2), Y(1,3), 'o', 'MarkerSize',15,'DisplayName','Firing start'), hold on
plot3(Y(end,1), Y(end,2), Y(end,3), 'ro', 'MarkerSize',15,'DisplayName','Firing end'), hold on

YY = [YY; Y];
TT = [TT; T];

figure()
delta_v_eff = parameters.Isp * 9.81 * log(parameters.M0./(parameters.M0 - YY(:,7)));
plot(TT, delta_v_eff), hold on, title('Delta v profile'), ylabel('Delta v [m/s]'), xlabel('t [s]')
yline(parameters.delta_v_req,'r'), hold on,
xline(T(1),'r'), hold on, xline(T(end),'r'), hold off

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


%% Thrust de-activated
parameters.t0sym = parameters.tmax;

T_orb = 2*pi*sqrt(YY_kep(end,1)^3/mu);
parameters.tmax = parameters.t0sym + T_orb/86400;
parameters.T = [0; 0; 0];            

[T, Y_arr] = cart_cont_thrust_model(YY(end,:), parameters);

TT = [TT; T];
figure1 = figure(1)
plot3(Y_arr(:,1),Y_arr(:,2),Y_arr(:,3),'b','DisplayName','Effective Capture Ellipse'), hold on
legend()

N_firings = 1;

annotation(figure1,'textbox',...
    [0.075 0.65 0.15 0.3],...
    'String', ...
    {'Desired Orbit Kep:', ...
    strcat('a: ', num2str(kep_cap(1)), 'km'), ...
    strcat('e: ', num2str(kep_cap(2))), ...
    strcat('i: ', num2str(rad2deg(kep_cap(3))), '°'), ...
    strcat('OM: ', num2str(rad2deg(kep_cap(4))), '°'), ...
    strcat('om: ', num2str(rad2deg(kep_cap(5))), '°'), ...
    strcat('theta: ', num2str(rad2deg(kep_cap(6))), '°'), ...
    ' ', ...
    'Effective Orbit Kep:', ...
    strcat('a: ', num2str(YY_kep(end,1)), 'km'), ...
    strcat('e: ', num2str(YY_kep(end,2))), ...
    strcat('i: ', num2str(rad2deg(YY_kep(end,3))), '°'), ...
    strcat('OM: ', num2str(rad2deg(YY_kep(end,4))), '°'), ...
    strcat('om: ', num2str(rad2deg(YY_kep(end,5))), '°'), ...
    strcat('theta: ', num2str(rad2deg(YY_kep(end,6))), '°'), ...
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