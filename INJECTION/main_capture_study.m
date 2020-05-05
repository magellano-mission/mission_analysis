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
Thrust0 = [-3200; 0; 0];                                   % thrust [N] (@TNH) (constant profile for the moment)
parameters.Isp = 280;                                      % specific impulse [s]
parameters.M0 = 5000;                                      % Total Mass of the s/c [kg]
parameters.c_r = 0.5;
parameters.Across_sun = 10;                                % Cross area related to the sun [m^2]
parameters.t0sym = date2mjd2000([2021, 1, 1, 0, 0, 0]);    
parameters.dt_p = parameters.t_BO/2;                       %time distance from pericenter at which firing occurs [s]
parameters.event = 0;
parameters.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12);   %,'Events', @event_cont_thrust);

%capture delta_v
kep_cap_desired = [15000 0.1 deg2rad(15) 0 0 0];
delta = deg2rad(45);
[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);
mu = astroConstants(14);

N_firings = 1;
[TT, YY, parout, dt_hyp, dv_req, theta_inf, kep_hyp_arr, kep_capture] = capture_plot(kep_cap_desired, delta,rM, mu, Thrust0, N_firings, parameters);
%the idea is to start the burn at (t_pericenter - t_BO/2) to emulate at
%most an impulsive maneuver at pericenter
%% non funziona
parameters.T = [0; 0; 0];                                  % thrust [N] (@TNH) (constant profile for the moment)

x0 = [parameters.dt_p];
% poi provo con il multiple firing
options = optimset('TolX',1e-10, 'TolFun', 1e-10);
x = fsolve(@(dt_p) single_burn_capture(Thrust0, parameters.Isp, dt_p, delta, rM, kep_cap_desired, mu, parameters), x0,options);

%%
[~, ~, ~, ~, ~, kep_hyp_arr, kep_capture] = capture_plot(kep_cap_desired, delta,rM, mu, x(1:3), N_firings, parameters);
%%

function [err] = single_burn_capture(Thrust, Isp, dt_p, delta, rM, kep_cap_desired, mu, parameters)
    parameters.T = [0; 0; 0];       % thrust [N] (@TNH) (constant profile for the moment)
    parameters.dt_p = dt_p;
    parameters.t0sym = date2mjd2000([2021, 1, 1, 0, 0, 0]);
    [~, ~, ~, ~, ~, ~, kep_cap_arr] = hyp2PO(kep_cap_desired, delta,rM, mu, Thrust , parameters);
    err = kep_cap_arr(1:2) - kep_cap_desired(1:2);
end