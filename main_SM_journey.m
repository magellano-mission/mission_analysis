%% INTERPLANETARY TRANSFER
%Propagation of the service module journey from PO (?) to Mars PO
%change the Thrust vector, M0, t_BO, Isp, (Mp0 is an output at the moment,
%myabe it is better to understand the order of magnitudes involved)
%delta_t_before_pericenter: trying to understand the impact of the firing
%timing

%these are all the parameters that we are using, changing among the code:
%t_BO
%thrust
%t0sym tmax

close, clear, clc

parameters.isInterp = 0;
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

%%%%%%  PHASES
% LEOP (lastly defined)
%% Interplanetary arc and TCM study

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

%%
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
R_M = R_E;
kep_E = zeros(length(T_interp),6);
kep_M = kep_E;
kep_interp = kep_E;

for jj = 1:length(T_interp)
    kep_E(jj,:) = uplanet(T_interp(jj)/86400,3);
    kep_M(jj,:) = uplanet(T_interp(jj)/86400,4);
    R_E(jj,:) = kep2car2(kep_E(jj,1:6), mu);
    R_M(jj,:) = kep2car2(kep_M(jj,1:6), mu);
    kep_interp(jj,:) = car2kep(Y_interp(jj,1:3), Y_interp(jj,4:6), mu);
end

%@cartesian r/f (unperturbed)
parameters.isPerturbed = 0;
[T_interp_unp, Y_interp_unp, parout_interp_unp] = cart_cont_thrust_model(X0, parameters);
for i=1:length(T_interp_unp)
    kep_interp_unp = car2kep(Y_interp_unp(i,1:3), Y_interp_unp(i,4:6), mu);
    kep_interp_unp(4:6) = rad2deg(kep_interp_unp(4:6));
end

figure()
hold all,
sgtitle('Keplerian parameters: Unperturbed vs Perturbed model')
txtkep{1} = 'a [km]'; txtkep{2} = 'e [-]'; txtkep{3} = 'i [deg]';
txtkep{4} = 'RAAN [deg]'; txtkep{5} = 'AP [deg]'; txtkep{6} = 'theta [deg]';
for i = 1:6
    subplot(3,3,i), plot(T_interp_unp, kep_interp_unp(:,i),'r:'), hold on
                    plot(T_interp, kep_interp(:,i),'b'),
                    title(txtkep(i))
end
legend('Unperturbed','Perturbed')
subplot(3,3,[7 8 9]), plot(T_interp,Y_interp(:,7),'b'); title('Expelled propellant mass')

% interplanetary arc plot 
figure_interp = figure;
I = imread('Sun.jpg'); RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];  RI.YWorldLimits = [-90 90]; 
rSun = almanac('Sun','Radius','kilometers','sphere');
[X, Y, Z] = ellipsoid(0, 0, 0, rSun, rSun, rSun, 100); % spheric centered Mars
planet = surf(X, Y, -Z,'Edgecolor', 'none'); hold on
set(planet,'FaceColor','texturemap','Cdata',I), axis equal

plot3(R_E(:,1), R_E(:,2), R_E(:,3),'b','DisplayName','Earth trajectory'), hold on
plot3(R_M(:,1), R_M(:,2), R_M(:,3),'r','DisplayName','Mars trajectory'), hold on

plot3(Y_interp_unp(1,:),Y_interp_unp(2,:),Y_interp_unp(3,:),'y:', 'DisplayName','Unperturbed'), hold on
plot3(Y_interp(1,:), Y_interp(2,:), Y_interp(3,:),'y', 'DisplayName','Perturbed'), hold on
plot3(Y_interp(1,1), Y_interp(2,1), Y_interp(3,1),'bo', 'DisplayName','Injection'), hold on
% set(gcf, 'color', 'k')
% set(gca, 'color', 'k','visible','off'), 
legend()
hold off 
% 
% figure()
% plot(T_interp, parameters.M0 - Y_interp(:,7)); title('s/c mass')

%% MOI and injection  
%Preliminary computations 
% Departure hyp
Ra_PO_dep = 315 + astroConstants(23);
Rp_PO_dep = 117 + astroConstants(23);

a_PO = (Ra_PO_dep + Rp_PO_dep)/2;
e_PO = (Ra_PO_dep - Rp_PO_dep)/(Ra_PO_dep + Rp_PO_dep);
p_PO = a_PO *(1-e_PO^2);

v_a_PO = sqrt( astroConstants(13)/p_PO )*(1-e_PO);
v_p_PO = sqrt( astroConstants(13)/p_PO )*(1+e_PO);
v_inf_dep = Dv1_best;
v_a_hyp = sqrt( ( v_inf_dep^2 ) + ( 2*astroConstants(13)/Ra_PO_dep ));
v_p_hyp = sqrt( ( v_inf_dep^2 ) + ( 2*astroConstants(13)/Rp_PO_dep ));
deltav_PO_a = v_a_hyp - v_a_PO;
deltav_PO_p = v_p_hyp - v_p_PO;

% Arrival hyp
Rp_PO_arrival = 11500:1000:15000;
v_PO_arr = sqrt( astroConstants(14)./Rp_PO_arrival );
v_inf_arrival = Dv2_best;
v_p_arr_hyp = sqrt( ( v_inf_arrival^2 ) + ( 2*astroConstants(14)./Rp_PO_arrival ));
deltav_PO_arr = v_p_arr_hyp - v_PO_arr;


