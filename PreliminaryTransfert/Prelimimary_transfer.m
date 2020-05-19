%% Figure Initialization
set(0,'DefaultFigureUnits', 'normalized');
set(0,'DefaultFigurePosition',[0 0 1 1]);
set(0,'DefaultTextFontSize',18);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultAxesXGrid','on')
set(0,'DefaultAxesYGrid','on')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

%% Preliminary Porkchop Plot

[~, mu_s] = uplanet(0, 1);

in_date_min = [2023 1 1 0 0 0];
fin_date_min = [2023 1 1 0 0 0];

[k_E, ~] = uplanet(date2mjd2000(in_date_min), 3);
[k_M, ~] = uplanet(date2mjd2000(in_date_min), 4);

E_P = 2*pi*sqrt(k_E(1)^3/mu_s);
M_P = 2*pi*sqrt(k_M(1)^3/mu_s);
n_per = 5;

% Synodic periods
EM_SP = (E_P*M_P) / abs(E_P-M_P);

[minDVI,Earth_time, Mars_time]=porkchop(in_date_min,fin_date_min,EM_SP,n_per);
colormap(MagellanoColorMap)

%% Deltav min for the IT, Delta V at the departure and arrival hyp
[k_E, ~] = uplanet(Earth_time, 3);
[r_E, v_E] = kep2car2([k_E(1), k_E(2), k_E(3), k_E(4), k_E(5), k_E(6)], mu_s);
        
[k_M, ~] = uplanet(Mars_time, 4);
[r_M, v_M] = kep2car2([k_M(1), k_M(2), k_M(3), k_M(4), k_M(5), k_M(6)], mu_s);
                
[~, ~ ,~, ~, VI, VF, ~, ~] = lambertMR(r_E, r_M, (Mars_time- Earth_time)*24*3600 , mu_s);
        
Dv1_best = norm(v_E - VI');
Dv2_best = norm(VF' - v_M);

% [~, Dv] = Earth_Mars_transfer_plot(Earth_time, Mars_time);

% Departure hyp
Ra_PO_dep = 315 + astroConstants(23);
% Rp_PO_dep = 117 + astroConstants(23);
Rp_PO_dep = (100:10:300) + astroConstants(23);

a_PO = (Ra_PO_dep + Rp_PO_dep)/2;
e_PO = (Ra_PO_dep - Rp_PO_dep)./(Ra_PO_dep + Rp_PO_dep);
p_PO = a_PO .*(1-e_PO.^2);

v_a_PO = sqrt( astroConstants(13)./p_PO ).*(1-e_PO);
v_p_PO = sqrt( astroConstants(13)./p_PO ).*(1+e_PO);
v_inf_dep = Dv1_best;
v_a_hyp = sqrt( ( v_inf_dep^2 ) + ( 2*astroConstants(13)/Ra_PO_dep ));
v_p_hyp = sqrt( ( v_inf_dep^2 ) + ( 2*astroConstants(13)./Rp_PO_dep ));
deltav_PO_a = v_a_hyp - v_a_PO;
deltav_PO_p = v_p_hyp - v_p_PO;

% Arrival hyp
Rp_PO_arrival = 11500:1000:15000;
v_PO_arr = sqrt( astroConstants(14)./Rp_PO_arrival );
v_inf_arrival = Dv2_best;
v_p_arr_hyp = sqrt( ( v_inf_arrival^2 ) + ( 2*astroConstants(14)./Rp_PO_arrival ));
deltav_PO_arr = v_p_arr_hyp - v_PO_arr;

figure()
plot(Rp_PO_dep - astroConstants(23), deltav_PO_p,'.-','DisplayName','Pericenter altitude'), hold on
plot(Rp_PO_dep - astroConstants(23), deltav_PO_p,'.-','DisplayName','Pericenter altitude')

%% Sun aspect angle during IT and around planetary orbit

% SAA e Earth occultation during the IT
SAA_tran = SAA_fun_tran(Earth_time, Mars_time);

% Earth aspect angle during IT
[EAA_] = Earth_visibility_angle(Earth_time, Mars_time);

% SAA in orbit around Mars
R_o = 11500;
mu_M = astroConstants(14);
[SAA_M] = SAA_fun(R_o,0,deg2rad(55),0,0,0,Mars_time,mu_M,4);

% SAA in orbit around Earth;
mu_E = astroConstants(13);
[SAA_E] = SAA_fun(a_PO(1),0,0,0,0,0,Earth_time,mu_E,3);

















