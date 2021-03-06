%% INTERPLANETARY TRANSFER AND TCM STUDY - LAUNCH 2
%Propagation of the service module journey from PO (?) to Mars PO
%change the Thrust vector, M0, t_BO, Isp, (Mp0 is an output at the moment,
%myabe it is better to understand the order of magnitudes involved)
%delta_t_before_pericenter: trying to understand the impact of the firing
%timing
close all, clear, clc
% Figure Initialization    
load('MagellanoColorMap.mat');
DefaultOrderColor = get(0, 'DefaultAxesColorOrder');
NewOrderColor = [0.9490    0.4745    0.3137
                 0.1020    0.6667    0.74120
                 155/255   155/255   155/255
                 DefaultOrderColor];  
             
set(0,'DefaultFigureColormap', MagellanoColorMap);
set(0, 'DefaultAxesColorOrder', NewOrderColor);
set(0,'DefaultLineLineWidth', 2)
set(0,'DefaultLineMarkerSize', 10)
set(0, 'DefaultFigureUnits', 'normalized');
set(0, 'DefaultFigurePosition', [0 0 1 1]);
set(0, 'DefaultTextFontSize', 18);
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

parameters.isPerturbed = 0;
parameters.isInterp = 1;
parameters.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
parameters.t_BO = 30*60;                                   % burnout time (chemical thrust)
parameters.T = [0; 0; 0];                                  % thrust [N] (@TNH) (constant profile for the moment)
Thrust0 = [-3200; 0; 0];                                   % thrust [N] (@TNH) (constant profile for the moment)
parameters.Isp = 280;                                      % specific impulse [s]
parameters.M0 = 5000;                                     % Total Mass of the s/c [kg]
parameters.c_r = 0.5;
parameters.Across_sun = 10;                                % Cross area related to the sun [m^2]
parameters.t0sym = date2mjd2000([2021, 1, 1, 0, 0, 0]);    
parameters.dt_p = parameters.t_BO/2;                       %time distance from pericenter at which firing occurs [s]
parameters.event = 0;
parameters.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12);   %,'Events', @event_cont_thrust);


%% Interplanetary arc - Lambert

% Preliminary Porkchop Plot

[~, mu_s] = uplanet(0, 1);

in_date_min = [2026 1 1 0 0 0];
fin_date_min = [2027 1 1 0 0 0];

[k_E, ~] = uplanet(date2mjd2000(in_date_min), 3);
[k_M, ~] = uplanet(date2mjd2000(in_date_min), 4);

E_P = 2*pi*sqrt(k_E(1)^3/mu_s);
M_P = 2*pi*sqrt(k_M(1)^3/mu_s);
n_per = 2;

% Synodic periods
EM_SP = (E_P*M_P) / abs(E_P-M_P);

[minDVI,Earth_time, Mars_time]=porkchop(in_date_min,fin_date_min,EM_SP,n_per);
Earth_date = mjd20002date(Earth_time);
Mars_date = mjd20002date(Mars_time);
annotation(gcf,'textarrow',[0.63125 0.46171875],...
    [0.454246913580247 0.54320987654321],...
    'String',{strcat('Departure: [', num2str(Earth_date(3)), '/',num2str(Earth_date(2)), '/', num2str(Earth_date(1)) ,']'), ...
    strcat('Arrival: [', num2str(Mars_date(3)), '/',num2str(Mars_date(2)), '/', num2str(Mars_date(1)) ,']'), ...
    strcat('Cost:', num2str(minDVI), 'km/s')});
load('MagellanoColorMap.mat')
colormap(MagellanoColorMap)
%%
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

% [~, Dv] = Earth_Mars_transfer_plot(Earth_time, Mars_time);

% %% Trans-Mars injection hyperbola 
% % https://www.nasa.gov/mission_pages/MRO/spacecraft/launch-cont.html
% 
% % First engine firing: The Centaur engine fires for the first time shortly 
% % after separation from the Stage I booster to boost the spacecraft into a PO of about 185 kilometers altitude
% % A "parking orbit" is the name used for the Earth orbit in which the spacecraft and Centaur coast between burns. 
% % % The parking orbit is needed to allow the Centaur to be in the right position relative to both Earth and Mars for each of its two burns. 
% % % The duration of the first burn is about nine-and-a-half minutes.
% % 
% 
% parameters.isPerturbed = 0;
% parameters.isInterp = 1;
% parameters.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
% parameters.T = [0; 0; 0];                                  % thrust [N] (@TNH) (constant profile for the moment)
% Thrust0 = [100000; 0; 0];                                   % thrust [N] (@TNH) (constant profile for the moment)
% parameters.Isp = 380;                                      % specific impulse [s]
% parameters.M0 = 14000;                                      % Total Mass of the s/c [kg]
% parameters.t0sym = (Earth_time);    
% parameters.event = 0;
% parameters.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12);   %,'Events', @event_cont_thrust);
% parameters.perc_tBO_before_p = 0.5;
% 
% 
% h_PO = 10000;
% rEarth = almanac('Earth','Radius','kilometers','sphere');
% r_PO = h_PO + rEarth;       %[km]
% muE = astroConstants(13); 
% v_inf_plus = v_E' - VI;   
% % rp_opt_parc = 2*muE/norm(v_inf_plus)^2;%TBC
% e_opt_parc = 0.1;
% a_opt_parc = r_PO/(1 - e_opt_parc);
% i_opt_parc = deg2rad(26.7);
% kep_parc = [a_opt_parc e_opt_parc i_opt_parc 0 0 0];%TBC
% 
% [ YM_NS, hyp_NS , kep_cap_NS] = PO2hyp(kep_parc, v_inf_plus, r_E, muE, parameters, Thrust0, 1, 'departure', 1);
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
R_M = R_E; R_J = R_M; R_V = R_M;
kep_E = zeros(length(T_interp),6);
kep_M = kep_E; kep_J = kep_M;  kep_V = kep_M;
kep_interp = kep_E;

for jj = 1:length(T_interp)
    kep_E(jj,:) = uplanet(T_interp(jj)/86400,3);
    kep_M(jj,:) = uplanet(T_interp(jj)/86400,4);
    kep_V(jj,:) = uplanet(T_interp(jj)/86400,2);
    kep_J(jj,:) = uplanet(T_interp(jj)/86400,5);
    R_E(jj,:) = kep2car2(kep_E(jj,1:6), mu);
    R_M(jj,:) = kep2car2(kep_M(jj,1:6), mu);
    R_V(jj,:) = kep2car2(kep_V(jj,1:6), mu);
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
planet = surf(X, Y, -Z,'Edgecolor', 'none','HandleVisibility','off'); hold on
set(planet,'FaceColor','texturemap','Cdata',I), axis equal

Marsplot = plot3(R_M(:,1), R_M(:,2), R_M(:,3),'DisplayName','Mars trajectory'); hold on
Earth_plot = plot3(R_E(:,1), R_E(:,2), R_E(:,3),'DisplayName','Earth trajectory'); hold on
% plot3(R_V(:,1), R_V(:,2), R_V(:,3),'c','DisplayName','Venus trajectory'), hold on
% plot3(R_J(:,1), R_J(:,2), R_J(:,3),'m','DisplayName','Jupiter trajectory'), hold on

% plot3(Y_interp_unp(:,1),Y_interp_unp(:,2),Y_interp_unp(:,3),'r', 'DisplayName','Unperturbed'), hold on
plot3(Y_interp(:,1), Y_interp(:,2), Y_interp(:,3),'Color', [150 150 150]/255, 'DisplayName','Perturbed'), hold on
xplot0 = plot3(Y_interp(1,1), Y_interp(1,2), Y_interp(1,3),'o', 'MarkerFaceColor', [0.1020    0.6667    0.74120], ...
    'MarkerEdgeColor', 'none', 'DisplayName','Injection'); hold on,
xplotend = plot3(Y_interp(end,1), Y_interp(end,2), Y_interp(end,3),'o', 'MarkerFaceColor', [0.9490    0.4745    0.3137], ...
    'MarkerEdgeColor', 'none', 'DisplayName','Capture'); hold on,


view([0 0 90])
% set(gcf, 'color', 'k')
% set(gca, 'color', 'k','visible','off'), 
legend('Color','white')

% hold off 
% 
% figure()
% plot(T_interp, parameters.M0 - Y_interp(:,7)); title('s/c mass')

% for i=1:10:length(T_interp)
%     rrr = [R_E(i,1:3); Y_interp(i,1:3)];
%     plot3(rrr(:,1), rrr(:,2), rrr(:,3), 'k','HandleVisibility','off'), hold on,
% end

% hold off 
% 
% figure()
% plot(T_interp, parameters.M0 - Y_interp(:,7)); title('s/c mass')
dist = zeros(length(T_interp),1);
ang = dist;
for i=1:length(dist)
    dist(i) = norm(Y_interp(i,1:3) - R_E(i,1:3));
    ang(i) = acosd(dot((Y_interp(i,1:3) - R_E(i,1:3))/norm((Y_interp(i,1:3) - R_E(i,1:3))), -R_E(i,1:3)/norm(R_E(i,1:3))));
end
% figure()
% % subplot(2,1,1), plot(T_interp, dist), xlabel('Mission Time [days]'), 
% % xlabel('Mission Time [days]'), ylabel('km'), title('Distance'), grid minor
% % subplot(2,1,2),
% SAA = plot((T_interp - T_interp(1))/86400, ang);
% % SAA.Color = [0.9490    0.4745    0.3137];
% SAA.Color = [0.1020    0.6667    0.74120];
% xlabel('Mission Time [days]'), ylabel('deg'), title('Angle Sun-Earth-s/c')
% grid minor, colormap(MagellanoColorMap)

%% TCM -> CAPTURE STUDY (2 SM)
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

parNS3 = parameters;
parRS = parameters;

dv_instantNS3 = ceil(0.95*length(T_interp));
TCM_before_MarsNS3 = (T_interp(end) - T_interp(dv_instantNS3))/86400 %days
TCM_dateNS3 = mjd20002date(T_interp(dv_instantNS3)/86400)
%capture delta_v
[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);
mu = astroConstants(14);

%second launch
kep_NS3 = [10500 0 deg2rad(25) 0 0 0]; 
kep_RS = [6400 0 deg2rad(0) 0 0 0];

%second launch
%NS
% [ YM_NS3, hyp_NS3, kep_cap_NS3] = ...
% capture_plot(kep_NS3, (VF' - v_M), rM, mu, Thrust0 , 1, parNS3);
[ YM_NS3, hyp_NS3, kep_cap_NS3] = PO2hyp(kep_NS3, (VF' - v_M), rM, mu, parNS3, [], 1, 'capture');

%definition of the last TCM
YSOI_NS3 = rM + YM_NS3(1, 1:3)';
[a_NS3, p_NS3 ,e_NS3, err_NS3, VI_NS3, VF_NS3, tspar_NS3, th_NS3] = lambertMR(Y_interp(dv_instantNS3,1:3)', YSOI_NS3, (Mars_time* 86400- T_interp(dv_instantNS3)) , mu_s);
dv_NS3 = VI_NS3 - Y_interp(dv_instantNS3, 4:6);
TCM_NS3 = norm(dv_NS3)*1000

[ ~, hyp_NS3, ~] = hyp2PO(kep_NS3, (VF_NS3' - v_M), rM, mu, parNS3, Thrust0, 1, 'capture');

parNS3.t0sym = T_interp(dv_instantNS3)/86400;
parNS3.tmax = Mars_time;
parNS3.isInterp = 1;
[~, Y_NS3, ~] = cart_cont_thrust_model([Y_interp(dv_instantNS3,1:3), VI_NS3, 0], parNS3);

% RS
dv_instantRS = ceil(0.95*length(T_interp));
TCM_before_MarsRS = (T_interp(end) - T_interp(dv_instantRS))/86400 %days
TCM_dateRS = mjd20002date(T_interp(dv_instantRS)/86400)
parRS.isInterp = 1;
% [ YM_RS, hyp_RS, kep_cap_RS] = ...
% capture_plot(kep_RS, (VF' - v_M), rM, mu, Thrust0 , 1, parRS);
[ YM_RS, hyp_RS, kep_cap_RS] = hyp2PO(kep_RS, (VF' - v_M), rM, mu, parRS, [], 1, 'capture');

%definition of the last TCM
YSOI_RS = rM + YM_RS(1, 1:3)';
[a_RS, p_RS ,e_RS, err_RS, VI_RS, VF_RS, tspar_RS, th_RS] = lambertMR(Y_interp(dv_instantRS,1:3)', YSOI_RS, (Mars_time*86400 - T_interp(dv_instantRS)) , mu_s);
dv_RS = VI_RS - Y_interp(dv_instantRS, 4:6);
[ ~, hyp_RS, ~] = hyp2PO(kep_RS, (VF_RS' - v_M), rM, mu, parRS, Thrust0, 1, 'capture');
TCM_RS = norm(dv_RS)*1000

parRS.t0sym = T_interp(dv_instantRS)/86400;
parRS.tmax = Mars_time;
parRS.isInterp = 1;
[~, Y_RS, ~] = cart_cont_thrust_model([Y_interp(dv_instantRS,1:3), VI_RS, 0], parRS);

% figure(3)
% plot3(Y_NS3(:,1), Y_NS3(:,2), Y_NS3(:,3), 'm', 'DisplayName','NS3 injection'); hold on
% plot3(Y_NS3(1,1), Y_NS3(1,2), Y_NS3(1,3), 'o', 'Color', 'm', ...
%     'DisplayName','TCM NS-3'); hold on
plot3(Y_RS(1,1), Y_RS(1,2), Y_RS(1,3),  'o', 'Color','k', 'DisplayName','Service Modules TCM'); hold on
plot3(Y_RS(:,1), Y_RS(:,2), Y_RS(:,3), 'k', 'DisplayName','NS+ECS and RS injection');hold on
view([0 0 90])
legend()
% annotation(figure_interp,'textbox', ...
%     [0.75 0.15 0.2 0.45], ...
%     'String',{'TCM:', ...
%     'NS3:', strcat('TCM:', num2str(1000*norm(dv_NS3)),' m/s', 'Capture: ',num2str(norm(hyp_NS3.dv_req)),' m/s', '(opt: ', num2str(norm(1000*hyp_NS3.dv_opt)),' m/s)') , ...
%     'RS:', strcat('TCM:', num2str(1000*norm(dv_RS)),' m/s', 'Capture: ',num2str(norm(hyp_RS.dv_req)),' m/s', '(opt: ', num2str(norm(1000*hyp_RS.dv_opt)),' m/s)')}, ...
%     'FitBoxToText','off');

%% TCM sensitivity study
parameters.isInterp = 0;
parameters.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
parameters.isPerturbed = 0;                                % perturbation switch 0 ,too slow :(
parameters.T = [0; 0; 0];                                  % thrust [N] (@TNH) (constant profile for the moment)
Thrust0 = [-10000; 0; 0];                                   % thrust [N] (@TNH) (constant profile for the moment)
parameters.t0sym = Mars_time;    
parameters.event = 0;
parameters.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12); 

%second launch
parNS3 = parameters;
parRS = parameters;
DV_NS3 = []; DV_RS = [];

perc_timing=0.01:0.01:0.99;


for tt = 1:length(perc_timing)
dv_instant = ceil(perc_timing(tt)*length(T_interp));

%capture delta_v
[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);
mu = astroConstants(14);

%NS3
% [ YM_NS3, ~, kep_cap_NS3] = ...
% capture_plot(kep_NS3, (VF' - v_M), rM, mu, Thrust0 , 1, parNS3);
[ YM_NS3, ~, kep_cap_NS3] = hyp2PO(kep_NS3, (VF' - v_M), rM, mu, parNS3, Thrust0, 1, 'capture');

%definition of the last TCM
YSOI_NS3 = rM + YM_NS3(1, 1:3)';
[~, ~ ,~, ~, VI_NS3, VF_NS3, ~, ~] = lambertMR(Y_interp(dv_instant,1:3)', YSOI_NS3, (Mars_time* 86400- T_interp(dv_instant)) , mu_s);
dv_NS3 = VI_NS3 - Y_interp(dv_instant, 4:6);

% [ ~, hyp_NS3, ~] = ...
% capture_plot(kep_NS3, (VF_NS3' - v_M), rM, mu, Thrust0 , 1, parNS3);
[ ~, hyp_NS3, ~] = hyp2PO(kep_NS3, (VF_NS3' - v_M), rM, mu, parNS3, Thrust0, 1, 'capture');

DV_NS3 = [DV_NS3; norm(dv_NS3), hyp_NS3.dv_req, hyp_NS3.dv_opt];

% RS
parRS.isInterp = 1;
% [ YM_RS, ~, ~] = ...
% capture_plot(kep_RS, (VF' - v_M), rM, mu, Thrust0 , 1, parRS);
[ YM_RS, ~, ~] = hyp2PO(kep_RS, (VF' - v_M), rM, mu, parRS, Thrust0, 1, 'capture');

%definition of the last TCM
YSOI_RS = rM + YM_RS(1, 1:3)';
[~, ~ ,~, ~, VI_RS, VF_RS, ~, ~] = lambertMR(Y_interp(dv_instant,1:3)', YSOI_RS, (Mars_time*86400 - T_interp(dv_instant)) , mu_s);
dv_RS = VI_RS - Y_interp(dv_instant, 4:6);

% [ ~, hyp_RS, ~] = ...
% capture_plot(kep_RS, (VF_RS' - v_M), rM, mu, Thrust0 , 1, parRS);
[ ~, hyp_RS, ~] = hyp2PO(kep_RS, (VF_RS' - v_M), rM, mu, parRS, Thrust0, 1, 'capture');

DV_RS = [DV_RS; norm(dv_RS), hyp_RS.dv_req, hyp_RS.dv_opt];
end

tof_interp = (T_interp(end) - T_interp(1)); % [s]
figure
% sgtitle('Launch 2: Single burn TCM and Capture Maneuver - Sensitivity to maneuver instant')
subplot(3,1,1), plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS3(:,1),'DisplayName', 'TCM NS3'), hold on
plot((1-perc_timing)*tof_interp/86400, 1000*DV_RS(:,1),'DisplayName', 'TCM RS'), hold on
xline((T_interp(end) - T_interp(dv_instantNS3))/86400, 'r', 'linewidth', 2, 'DisplayName', 'Selected TCM')
xlabel('Time before arriving on Mars [days]'), ylabel('\Deltav [m/s]'), grid minor, title('TCM maneuver')
legend()

subplot(3,1,2), plot((1-perc_timing)*tof_interp/86400, DV_NS3(:,2),'DisplayName', 'dv NS3'), hold on
plot((1-perc_timing)*tof_interp/86400, DV_RS(:,2),'DisplayName', 'dv RS'), hold on
xline((T_interp(end) - T_interp(dv_instantNS3))/86400, 'r', 'linewidth', 2, 'DisplayName', 'Selected TCM')
% plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS3(:,3), 's','DisplayName', 'opt dv NS3'), hold on
% plot((1-perc_timing)*tof_interp/86400, 1000*DV_RS(:,3),'s','DisplayName', 'opt dv RS'), hold on
xlabel('Time before arriving on Mars [days]'), ylabel('\Deltav [m/s]'), grid minor, title('Capture \Deltav')
legend()

subplot(3,1,3), plot((1-perc_timing)*tof_interp/86400, 1000*DV_NS3(:,1) + DV_NS3(:,2),'DisplayName', 'dv NS3'), hold on
plot((1-perc_timing)*tof_interp/86400, 1000*DV_RS(:,1) + DV_RS(:,2),'DisplayName', 'dv RS'), hold on
xline((T_interp(end) - T_interp(dv_instantNS3))/86400, 'r', 'linewidth', 2, 'DisplayName', 'Selected TCM')
xlabel('Time before arriving on Mars [days]'), ylabel('\Deltav [m/s]'), grid minor, title('Total \Deltav')
legend()


