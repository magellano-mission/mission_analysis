%% INTERPLANETARY TRANSFER AND TCM STUDY
%Propagation of the service module journey from PO (?) to Mars PO
%change the Thrust vector, M0, t_BO, Isp, (Mp0 is an output at the moment,
%myabe it is better to understand the order of magnitudes involved)
%delta_t_before_pericenter: trying to understand the impact of the firing
%timing

%these are all the parameters that we are using, changing among the code:
%t_BO
%thrust
%t0sym tmax

close all, clear, clc

parameters.isInterp = 1;
parameters.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
parameters.t_BO = 30*60;                                   % burnout time (chemical thrust)
parameters.T = [0; 0; 0];                                  % thrust [N] (@TNH) (constant profile for the moment)
Thrust0 = [-3200; 0; 0];                                   % thrust [N] (@TNH) (constant profile for the moment)
parameters.Isp = 280;                                      % specific impulse [s]
parameters.M0 = 5000;                                      % Total Mass of the s/c [kg]
parameters.c_r = 0.5;
parameters.Across_sun = 10;                                % Cross area related to the sun [m^2]
parameters.t0sym = date2mjd2000([2021, 1, 1, 0, 0, 0]);    
parameters.dt_p = parameters.t_BO/2;                       %time distance from pericenter at which firing occurs [s]
parameters.event = 0;
parameters.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12);   %,'Events', @event_cont_thrust);


%% Interplanetary arc - Lambert

% Preliminary Porkchop Plot

[~, mu_s] = uplanet(0, 1);

in_date_min = [2022 1 1 0 0 0];
fin_date_min = [2022 1 1 0 0 0];

[k_E, ~] = uplanet(date2mjd2000(in_date_min), 3);
[k_M, ~] = uplanet(date2mjd2000(in_date_min), 4);

E_P = 2*pi*sqrt(k_E(1)^3/mu_s);
M_P = 2*pi*sqrt(k_M(1)^3/mu_s);
n_per = 5;

% Synodic periods
EM_SP = (E_P*M_P) / abs(E_P-M_P);

[minDVI,Earth_time, Mars_time]=porkchop(in_date_min,fin_date_min,EM_SP,n_per);

parameters.t0sym = date2mjd2000([2021, 1, 1, 0, 0, 0]);    
parameters.tmax = date2mjd2000([2021, 1, 1, 0, 0, 0]);    

% Deltav min for the IT, Delta V at the departure and arrival hyp
[k_E, ~] = uplanet(Earth_time, 3);
[r_E, v_E] = kep2car2([k_E(1), k_E(2), k_E(3), k_E(4), k_E(5), k_E(6)], mu_s);
        
[k_M, ~] = uplanet(Mars_time, 4);
[r_M, v_M] = kep2car2([k_M(1), k_M(2), k_M(3), k_M(4), k_M(5), k_M(6)], mu_s);
                
[~, ~ ,~, ~, VI, VF, ~, ~] = lambertMR(r_E, r_M, (Mars_time- Earth_time)*24*3600 , mu_s);
        
Dv1_best = norm(v_E - VI');
Dv2_best = norm(VF' - v_M);

[~, Dv] = Earth_Mars_transfer_plot(Earth_time, Mars_time);

%% Trans-Mars injection hyperbola 

muE = astroConstants(13); 
v_inf_plus = v_E - VI;
rp_opt_parc = 2*muE/norm(v_inf_plus)^2;%TBC
e_opt_parc = 0;
a_opt_parc = rp_opt_parc/(1 - e_opt_parc);
i_opt_parc = deg2rad(26.7);
kep_parc = [a_opt_parc e_opt_parc i_opt_parc 0 0 0];%TBC

% [ YM_NS, hyp_NS , kep_cap_NS] = PO2hyp(kep_NS, (VF' - v_M), rM, mu, parNS, [], 1, 'capture');

%% Propagation of perturbed model
%perturbed model (and later TCM maneuvers definition)
%initialization

parameters.isInterp = 1;
parameters.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
parameters.isPerturbed = 1;                                % perturbation switch 0 ,too slow :(
parameters.t_BO = 30*60;                                   % burnout time (chemical thrust)
parameters.T = [0; 0; 0];                                  % thrust [N] (@TNH) (constant profile for the moment)
Thrust0 = [0; 0; 0]; %engine for TCM tbd                   % thrust [N] (@TNH) (constant profile for the moment)
parameters.Isp = 280;                                      % specific impulse [s]
parameters.M0 = 5000;                                      % Total Mass of the s/c [kg]

parameters.t0sym = Earth_time;    
parameters.tmax = Mars_time;

parameters.event = 0;
parameters.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12);   %,'Events', @event_cont_thrust);

mu = astroConstants(4);
        
X0 = zeros(7,1);
X0(1:3) = r_E; X0(4:6) = VI;

%@cartesian r/f
[T_interp, Y_interp, parout_interp] = cart_cont_thrust_model(X0, parameters);

%plot
R_E = zeros(length(T_interp),3);
R_M = R_E; R_J = R_M;
kep_E = zeros(length(T_interp),6);
kep_M = kep_E; kep_J = kep_M;
kep_interp = kep_E;

for jj = 1:length(T_interp)
    kep_E(jj,:) = uplanet(T_interp(jj)/86400,3);
    kep_M(jj,:) = uplanet(T_interp(jj)/86400,4);
    kep_J(jj,:) = uplanet(T_interp(jj)/86400,5);
    R_E(jj,:) = kep2car2(kep_E(jj,1:6), mu);
    R_M(jj,:) = kep2car2(kep_M(jj,1:6), mu);
    R_J(jj,:) = kep2car2(kep_J(jj,1:6), mu);
    kep_interp(jj,:) = car2kep(Y_interp(jj,1:3), Y_interp(jj,4:6), mu);
end

%@cartesian r/f (unperturbed)
parameters.isPerturbed = 0;
[T_interp_unp, Y_interp_unp, parout_interp_unp] = cart_cont_thrust_model(X0, parameters);
kep_interp_unp = zeros(length(T_interp_unp), 6);
for i=1:length(T_interp_unp)
    kep_interp_unp(i,:) = car2kep(Y_interp_unp(i,1:3), Y_interp_unp(i,4:6), mu);
    kep_interp_unp(i,4:6) = rad2deg(kep_interp_unp(i,4:6));
end

figure()
hold all,
sgtitle('Keplerian parameters: Unperturbed vs Perturbed model')
txtkep{1} = 'a [km]'; txtkep{2} = 'e [-]'; txtkep{3} = 'i [deg]';
txtkep{4} = 'RAAN [deg]'; txtkep{5} = 'AP [deg]'; txtkep{6} = 'theta [deg]';

for i = 1:6
%     kep_interp_unp1 = interp1(T_interp_unp, kep_interp_unp(:,i)', T_interp ,'nearest');
    subplot(3,3,i), plot(T_interp_unp, kep_interp_unp(:,i),'r:'), hold on
    plot(T_interp, kep_interp(:,i) ,'b'), hold on
    title(txtkep(i))
end
legend('Unperturbed','Perturbed')
subplot(3,3,[7 8 9]), plot(T_interp,Y_interp(:,7),'b'); title('Expelled propellant mass')

% interplanetary arc plot 
figure_interp = figure;
I = imread('Sun.jpg'); RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];  RI.YWorldLimits = [-90 90]; 
rSun = 20*almanac('Sun','Radius','kilometers','sphere');
[X, Y, Z] = ellipsoid(0, 0, 0, rSun, rSun, rSun, 100); % spheric centered Mars
planet = surf(X, Y, -Z,'Edgecolor', 'none','DisplayName', 'Sun'); hold on
set(planet,'FaceColor','texturemap','Cdata',I), axis equal

plot3(R_E(:,1), R_E(:,2), R_E(:,3),'b','DisplayName','Earth trajectory'), hold on
plot3(R_M(:,1), R_M(:,2), R_M(:,3),'r','DisplayName','Mars trajectory'), hold on
plot3(R_J(:,1), R_J(:,2), R_J(:,3),'r','DisplayName','Jupiter trajectory'), hold on

plot3(Y_interp_unp(:,1),Y_interp_unp(:,2),Y_interp_unp(:,3),'r:', 'DisplayName','Unperturbed'), hold on
plot3(Y_interp(:,1), Y_interp(:,2), Y_interp(:,3),'g', 'DisplayName','Perturbed'), hold on
plot3(Y_interp(1,1), Y_interp(1,2), Y_interp(1,3),'go', 'DisplayName','Injection'), hold on,
plot3(Y_interp(end,1), Y_interp(end,2), Y_interp(end,3),'ro', 'DisplayName','Capture'), hold on,

view([0 0 90])
% set(gcf, 'color', 'k')
% set(gca, 'color', 'k','visible','off'), 
legend('Color','white')

% hold off 
% 
% figure()
% plot(T_interp, parameters.M0 - Y_interp(:,7)); title('s/c mass')

%% TCM -> CAPTURE STUDY (4 SM)
%these are all the parameters that we are using, changing among the code:
%t_BO
%thrust
%t0sym tmax

parameters.isInterp = 0;
parameters.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
parameters.isPerturbed = 0;                                % perturbation switch 0 ,too slow :(
parameters.T = [0; 0; 0];                                  % thrust [N] (@TNH) (constant profile for the moment)
Thrust0 = [-10000; 0; 0];                                   % thrust [N] (@TNH) (constant profile for the moment)
parameters.t0sym = Mars_time;    
parameters.event = 0;
parameters.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12); 

parNS = parameters;
parNS2 = parameters;
parNS3 = parameters;
parRS = parameters;

dv_instant = ceil(0.8*length(T_interp));

%capture delta_v
[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);
mu = astroConstants(14);

%first launch
kep_NS = [9850 0 deg2rad(25) 0 0 0];
kep_NS2 = [8100 0 deg2rad(25) 0 0 0];
%second launch
kep_NS3 = [10500 0 deg2rad(25) 0 0 0]; 
kep_RS = [6400 0 deg2rad(0) 0 0 0];

%first launch
%NS
% [ YM_NS, hyp_NS, kep_cap_NS] = ...
% capture_plot(kep_NS, (VF' - v_M), rM, mu, Thrust0 , 1, parNS);
[ YM_NS, hyp_NS , kep_cap_NS] = PO2hyp(kep_NS, (VF' - v_M), rM, mu, parNS, [], 1, 'capture');

%definition of the last TCM
YSOI_NS = rM + YM_NS(1:3, 1);
[a_NS, p_NS ,e_NS, err_NS, VI_NS, VF_NS, tspar_NS, th_NS] = lambertMR(Y_interp(dv_instant,1:3)', YSOI_NS, (Mars_time* 86400- T_interp(dv_instant)) , mu_s);
dv_NS = VI_NS - Y_interp(dv_instant, 4:6);

parNS.t0sym = T_interp(dv_instant)/86400;
parNS.tmax = Mars_time;
parNS.isInterp = 1;
[~, Y_NS, ~] = cart_cont_thrust_model([Y_interp(dv_instant,1:3), VI_NS, 0], parNS);

%NS2
parNS2.isInterp = 1;
% [ YM_NS2, hyp_NS2, kep_cap_NS2] = ...
% capture_plot(kep_NS2, (VF' - v_M), rM, mu, Thrust0 , 1, parNS2);
[ YM_NS2, hyp_NS2, kep_cap_NS2] = PO2hyp(kep_NS2, (VF' - v_M), rM, mu, parNS2, [], 1, 'capture');

%definition of the last TCM
YSOI_NS2 = rM + YM_NS2(1:3, 1);
[a_NS2, p_NS2 ,e_NS2, err_NS2, VI_NS2, VF_NS2, tspar_NS2, th_NS2] = lambertMR(Y_interp(dv_instant,1:3)', YSOI_NS2, (Mars_time*86400 - T_interp(dv_instant)) , mu_s);
dv_NS2 = VI_NS2 - Y_interp(dv_instant, 4:6);

parNS2.t0sym = T_interp(dv_instant)/86400;
parNS2.tmax = Mars_time;
parNS2.isInterp = 1;
[~, Y_NS2, ~] = cart_cont_thrust_model([Y_interp(dv_instant,1:3), VI_NS2, 0], parNS2);

%second launch
%NS
% [ YM_NS3, hyp_NS3, kep_cap_NS3] = ...
% capture_plot(kep_NS3, (VF' - v_M), rM, mu, Thrust0 , 1, parNS3);
[ YM_NS3, hyp_NS3, kep_cap_NS3] = PO2hyp(kep_NS3, (VF' - v_M), rM, mu, parNS3, [], 1, 'capture');

%definition of the last TCM
YSOI_NS3 = rM + YM_NS3(1:3, 1);
[a_NS3, p_NS3 ,e_NS3, err_NS3, VI_NS3, VF_NS3, tspar_NS3, th_NS3] = lambertMR(Y_interp(dv_instant,1:3)', YSOI_NS3, (Mars_time* 86400- T_interp(dv_instant)) , mu_s);
dv_NS3 = VI_NS3 - Y_interp(dv_instant, 4:6);

parNS3.t0sym = T_interp(dv_instant)/86400;
parNS3.tmax = Mars_time;
parNS3.isInterp = 1;
[~, Y_NS3, ~] = cart_cont_thrust_model([Y_interp(dv_instant,1:3), VI_NS3, 0], parNS3);

% RS
parRS.isInterp = 1;
% [ YM_RS, hyp_RS, kep_cap_RS] = ...
% capture_plot(kep_RS, (VF' - v_M), rM, mu, Thrust0 , 1, parRS);
[ YM_RS, hyp_RS, kep_cap_RS] = PO2hyp(kep_RS, (VF' - v_M), rM, mu, parRS, [], 1, 'capture');

%definition of the last TCM
YSOI_RS = rM + YM_RS(1:3, 1);
[a_RS, p_RS ,e_RS, err_RS, VI_RS, VF_RS, tspar_RS, th_RS] = lambertMR(Y_interp(dv_instant,1:3)', YSOI_RS, (Mars_time*86400 - T_interp(dv_instant)) , mu_s);
dv_RS = VI_RS - Y_interp(dv_instant, 4:6);

parRS.t0sym = T_interp(dv_instant)/86400;
parRS.tmax = Mars_time;
parRS.isInterp = 1;
[~, Y_RS, ~] = cart_cont_thrust_model([Y_interp(dv_instant,1:3), VI_RS, 0], parRS);

figure(3)
plot3(Y_NS(:,1), Y_NS(:,2), Y_NS(:,3), 'k', 'DisplayName','NS injection'); hold on
plot3(Y_NS(1,1), Y_NS(1,2), Y_NS(1,3), 'ko', 'DisplayName','NS TCM');hold on
plot3(Y_NS2(:,1), Y_NS2(:,2), Y_NS2(:,3), 'k', 'DisplayName','NS injection');hold on
plot3(Y_NS2(1,1), Y_NS2(1,2), Y_NS2(1,3), 'ko', 'DisplayName','NS TCM');hold on
plot3(Y_NS3(:,1), Y_NS3(:,2), Y_NS3(:,3), 'k', 'DisplayName','NS3 injection'); hold on
plot3(Y_NS3(1,1), Y_NS3(1,2), Y_NS3(1,3), 'ko', 'DisplayName','NS3 TCM');hold on
plot3(Y_RS(:,1), Y_RS(:,2), Y_RS(:,3), 'k', 'DisplayName','NS injection');hold on
plot3(Y_RS(1,1), Y_RS(1,2), Y_RS(1,3), 'ko', 'DisplayName','NS TCM');hold on


annotation(figure_interp,'textbox', ...
    [0.75 0.15 0.2 0.45], ...
    'String',{'TCM:', ...
    'NS:', strcat('TCM:', num2str(1000*norm(dv_NS)),' m/s', 'Capture: ',num2str(norm(hyp_NS.dv_req)),' m/s', '(opt: ', num2str(norm(1000*hyp_NS.dv_opt)),' m/s)') , ...
    'NS2:', strcat('TCM:', num2str(1000*norm(dv_NS2)),' m/s', 'Capture: ',num2str(norm(hyp_NS2.dv_req)),' m/s', '(opt: ', num2str(norm(1000*hyp_NS2.dv_opt)),' m/s)'), ...
    'RS:', strcat('TCM:', num2str(1000*norm(dv_RS)),' m/s', 'Capture: ',num2str(norm(hyp_RS.dv_req)),' m/s', '(opt: ', num2str(norm(1000*hyp_RS.dv_opt)),' m/s)')}, ...
    'FitBoxToText','off');

%% TCM sensitivity study
parameters.isInterp = 0;
parameters.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
parameters.isPerturbed = 0;                                % perturbation switch 0 ,too slow :(
parameters.T = [0; 0; 0];                                  % thrust [N] (@TNH) (constant profile for the moment)
Thrust0 = [-10000; 0; 0];                                   % thrust [N] (@TNH) (constant profile for the moment)
parameters.t0sym = Mars_time;    
parameters.event = 0;
parameters.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12); 

%first launch
parNS = parameters;
parNS2 = parameters;
%second launch
parNS3 = parameters;
parRS = parameters;
DV_NS = []; DV_NS2 = []; 
DV_NS3 = []; DV_RS = [];

perc_timing=0.01:0.01:0.99;

%first launch
for tt = 1:length(perc_timing)
dv_instant = ceil(perc_timing(tt)*length(T_interp));

%capture delta_v
[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);
mu = astroConstants(14);

%NS
% [ YM_NS, ~, kep_cap_NS] = ...
% capture_plot(kep_NS, (VF' - v_M), rM, mu, Thrust0 , 1, parNS);
[ YM_NS, ~, kep_cap_NS] = PO2hyp(kep_NS, (VF' - v_M), rM, mu, parNS, Thrust0, 1, 'capture');

%definition of the last TCM
YSOI_NS = rM + YM_NS(1:3, 1);
[~, ~ ,~, ~, VI_NS, VF_NS, ~, ~] = lambertMR(Y_interp(dv_instant,1:3)', YSOI_NS, (Mars_time* 86400- T_interp(dv_instant)) , mu_s);
dv_NS = VI_NS - Y_interp(dv_instant, 4:6);

% [ ~, hyp_NS, ~] = ...
% capture_plot(kep_NS, (VF_NS' - v_M), rM, mu, Thrust0 , 1, parNS);
[ ~, hyp_NS, ~] = PO2hyp(kep_NS, (VF_NS' - v_M), rM, mu, parNS, Thrust0, 1, 'capture');

DV_NS = [DV_NS; norm(dv_NS), hyp_NS.dv_req, hyp_NS.dv_opt];

%NS2
parNS2.isInterp = 1;
% [ YM_NS2, ~, ~] = ...
% capture_plot(kep_NS2, (VF' - v_M), rM, mu, Thrust0 , 1, parNS2);
[ YM_NS2, ~, ~] = PO2hyp(kep_NS2, (VF' - v_M), rM, mu, parNS2, Thrust0, 1, 'capture');

%definition of the last TCM
YSOI_NS2 = rM + YM_NS2(1:3, 1);
[~, ~ ,~, ~, VI_NS2, VF_NS2, ~, ~] = lambertMR(Y_interp(dv_instant,1:3)', YSOI_NS2, (Mars_time*86400 - T_interp(dv_instant)) , mu_s);
dv_NS2 = VI_NS2 - Y_interp(dv_instant, 4:6);

% [ ~, hyp_NS2, ~] = ...
% capture_plot(kep_NS2, (VF_NS2' - v_M), rM, mu, Thrust0 , 1, parNS2);
[ ~, hyp_NS2, ~] = PO2hyp(kep_NS2, (VF_NS2' - v_M), rM, mu, parNS2, Thrust0, 1, 'capture');

DV_NS2 = [DV_NS2; norm(dv_NS2), hyp_NS2.dv_req, hyp_NS2.dv_opt];

end

tof_interp = (T_interp(end) - T_interp(1)); % [s]

figure
sgtitle('Launch 1: Single burn TCM and Capture Maneuver - Sensitivity to maneuver instant')
subplot(3,1,1), plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS(:,1),'DisplayName', 'TCM NS'), hold on
plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS2(:,1),'DisplayName', 'TCM NS2'), hold on
xlabel('Time before arriving on Mars [days]'), ylabel('delta v [m/s]'), grid minor, title('TCM maneuver')
legend()

subplot(3,1,2), plot((1-perc_timing)*tof_interp/86400, DV_NS(:,2),'DisplayName', 'dv NS'), hold on
plot((1-perc_timing)*tof_interp/86400, DV_NS2(:,2),'DisplayName', 'dv NS2'), hold on
plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS(:,3), 's','DisplayName', 'opt dv NS'), hold on
plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS2(:,3),'s','DisplayName', 'opt dv NS2'), hold on
xlabel('Time before arriving on Mars [days]'), ylabel('delta v [m/s]'), grid minor, title('Capture delta v')
legend()

subplot(3,1,3), plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS(:,1) + DV_NS(:,2),'DisplayName', 'dv NS'), hold on
plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS2(:,1) + DV_NS2(:,2),'DisplayName', 'dv NS2'), hold on
xlabel('Time before arriving on Mars [days]'), ylabel('delta v [m/s]'), grid minor, title('Total delta v')
legend()

%% second launch

for tt = 1:length(perc_timing)
dv_instant = ceil(perc_timing(tt)*length(T_interp));

%capture delta_v
[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);
mu = astroConstants(14);

%NS3
% [ YM_NS3, ~, kep_cap_NS3] = ...
% capture_plot(kep_NS3, (VF' - v_M), rM, mu, Thrust0 , 1, parNS3);
[ YM_NS3, ~, kep_cap_NS3] = PO2hyp(kep_NS3, (VF' - v_M), rM, mu, parNS3, Thrust0, 1, 'capture');

%definition of the last TCM
YSOI_NS3 = rM + YM_NS3(1:3, 1);
[~, ~ ,~, ~, VI_NS3, VF_NS3, ~, ~] = lambertMR(Y_interp(dv_instant,1:3)', YSOI_NS3, (Mars_time* 86400- T_interp(dv_instant)) , mu_s);
dv_NS3 = VI_NS3 - Y_interp(dv_instant, 4:6);

% [ ~, hyp_NS3, ~] = ...
% capture_plot(kep_NS3, (VF_NS3' - v_M), rM, mu, Thrust0 , 1, parNS3);
[ ~, hyp_NS3, ~] = PO2hyp(kep_NS3, (VF_NS3' - v_M), rM, mu, parNS3, Thrust0, 1, 'capture');

DV_NS3 = [DV_NS3; norm(dv_NS3), hyp_NS3.dv_req, hyp_NS3.dv_opt];

% RS
parRS.isInterp = 1;
% [ YM_RS, ~, ~] = ...
% capture_plot(kep_RS, (VF' - v_M), rM, mu, Thrust0 , 1, parRS);
[ YM_RS, ~, ~] = PO2hyp(kep_RS, (VF' - v_M), rM, mu, parRS, Thrust0, 1, 'capture');

%definition of the last TCM
YSOI_RS = rM + YM_RS(1:3, 1);
[~, ~ ,~, ~, VI_RS, VF_RS, ~, ~] = lambertMR(Y_interp(dv_instant,1:3)', YSOI_RS, (Mars_time*86400 - T_interp(dv_instant)) , mu_s);
dv_RS = VI_RS - Y_interp(dv_instant, 4:6);

% [ ~, hyp_RS, ~] = ...
% capture_plot(kep_RS, (VF_RS' - v_M), rM, mu, Thrust0 , 1, parRS);
[ ~, hyp_RS, ~] = PO2hyp(kep_RS, (VF_RS' - v_M), rM, mu, parRS, Thrust0, 1, 'capture');

DV_RS = [DV_RS; norm(dv_RS), hyp_RS.dv_req, hyp_RS.dv_opt];
end

tof_interp = (T_interp(end) - T_interp(1)); % [s]

figure
sgtitle('Launch 2: Single burn TCM and Capture Maneuver - Sensitivity to maneuver instant')
subplot(3,1,1), plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS3(:,1),'DisplayName', 'TCM NS3'), hold on
plot((1-perc_timing)*tof_interp/86400, 1000*DV_RS(:,1),'DisplayName', 'TCM RS'), hold on
xlabel('Time before arriving on Mars [days]'), ylabel('delta v [m/s]'), grid minor, title('TCM maneuver')
legend()

subplot(3,1,2), plot((1-perc_timing)*tof_interp/86400, DV_NS3(:,2),'DisplayName', 'dv NS3'), hold on
plot((1-perc_timing)*tof_interp/86400, DV_RS(:,2),'DisplayName', 'dv RS'), hold on
plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS3(:,3), 's','DisplayName', 'opt dv NS3'), hold on
plot((1-perc_timing)*tof_interp/86400, 1000*DV_RS(:,3),'s','DisplayName', 'opt dv RS'), hold on
xlabel('Time before arriving on Mars [days]'), ylabel('delta v [m/s]'), grid minor, title('Capture delta v')
legend()

subplot(3,1,3), plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS3(:,1) + DV_NS3(:,2),'DisplayName', 'dv NS3'), hold on
plot((1-perc_timing)*tof_interp/86400, 1000*DV_RS(:,1) + DV_RS(:,2),'DisplayName', 'dv RS'), hold on
xlabel('Time before arriving on Mars [days]'), ylabel('delta v [m/s]'), grid minor, title('Total delta v')
legend()

%% MOI injection
% B-plane
% figure('Name','B-plane')
% mu = astroConstants(14);
% BBB = [];
% for rp_desired = 8500:1000:15000
% for i_cap = 1:360
%     i_cap = deg2rad(i_cap);
% 
%     [e_cap, dv_min, delta_opt] = hyp_opt(rp_desired, (VF' - v_M), mu);
% 
%     kep_cap_desired = [rp_desired/(1-e_cap) e_cap i_cap 0 0 0];
% 
%     muS = astroConstants(4);
%     %Discretization of hyperbola
%     n_iter = 1000;
% 
%     delta = asin(1/(1 + rp_desired*norm(VF'-v_M)/mu ));
% 
%     %Computation of the hyperbola
%     e_m = 1/sin(delta/2);
%     h_m = sqrt(mu*rp_desired*(1+e_m));
%     a_m =((h_m^2)/mu)/(e_m^2-1);
%     theta_inf = acos(1/e_m);
%     nrM = norm(r_M); 
%     %SOI definition
%     rSOI = nrM * (mu/muS)^(2/5);
% 
%     correc_m_f = @(alpha) correc_SOI_f(rSOI,alpha,-a_m,e_m,theta_inf,mu);
%     correc_m = fminsearch(correc_m_f,0);
%     ltheta_m = linspace(0,theta_inf+correc_m,n_iter);
%     lt_m = zeros(1,n_iter);
%     %TOF computation
%     for i = 1:n_iter
%         %eccentric anomalies
%         f_m = 2*atanh(sqrt((e_m-1)/(e_m+1))*tan(ltheta_m(i)/2));
%         m_m = e_m*sinh(f_m) - f_m;
%         lt_m(i) = m_m*(h_m^3/mu^2)/((e_m^2-1)^(3/2));
%     end
%     %Computation of TOF
%     dt_hyp = abs(lt_m(1) - lt_m(n_iter)); % [s]
% 
%     %Hyperbolas plot
%     l_pos_Vect_m = zeros(3,n_iter);
% 
%     kep_hyp_arr = kep_cap_desired;
%     kep_hyp_arr(1) = -a_m;
%     kep_hyp_arr(2) = e_m;
% 
%     [~, v_ell] = kep2car2(kep_cap_desired, mu);
%     [~, v_hyp] = kep2car2(kep_hyp_arr, mu);
%     dv_req = v_ell - v_hyp; %just a check
% 
%     for i = 1:n_iter
%         kep_hyp_arr(6) = ltheta_m(i);
%         [r_m,~]= kep2car2(kep_hyp_arr,mu);
%         l_pos_Vect_m(:,i) = r_m;
%     end
% 
%     X0_hyp = zeros(6,1);
%     [rr_inf, vv_inf] = kep2car2(kep_hyp_arr, mu);
%     X0_hyp(1:3) = rr_inf;
%     X0_hyp(4:6) = -vv_inf;   
% 
%     [BB, B, theta, STR] = B_plane(kep_hyp_arr, X0_hyp(1:3), X0_hyp(4:6), mu);
% 
%     BBB = [BBB, BB];
%     end
% 
%     quiver3(0, 0, 0, STR(1,1), STR(2,1), STR(3,1), norm(BB)), hold on
%     quiver3(0, 0, 0, STR(1,2), STR(2,2), STR(3,2), norm(BB)), hold on
%     quiver3(0, 0, 0, STR(1,3), STR(2,3), STR(3,3), norm(BB)), hold on
% 
%     plot3(BBB(1,:), BBB(2,:), BBB(3,:)), hold on
% 
%     I = imread('Mars.jpg'); RI = imref2d(size(I));
%     RI.XWorldLimits = [-180 180];  RI.YWorldLimits = [-90 90]; 
%     rMars = almanac('Mars','Radius','kilometers','sphere');
%     [X, Y, Z] = ellipsoid(0, 0, 0, rMars, rMars, rMars, 100); % spheric centered Mars
%     planet = surf(X, Y, -Z,'Edgecolor', 'none','DisplayName', 'Sun'); hold on
%     set(planet,'FaceColor','texturemap','Cdata',I), axis equal
% 
%     legend('S', 'T', 'R')
%     % Function used to get a correction on the theta for the hyperbola arc to
%     % reach the Sphere of Influence
% end
function val = correc_SOI_f(rSOI,alpha,a,e,theta,muM)
    [r,~]= kep2car2([a,e,0,0,0,theta+alpha],muM);
    val = (norm(r)-rSOI)^2;
end
