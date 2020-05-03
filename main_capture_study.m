%% CAPTURE STUDY

close all, clear all, clc

parameters.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
parameters.t_BO = 3000000;                                 % burnout time (chemical thrust)
parameters.T = [-10; 0; 0];                                % thrust [N] (constant for the moment)
parameters.Is = 200;                                       % specific impulse [s]
parameters.M0 = 5000;                                      % Total Mass of the s/c [kg]
parameters.c_r = 0.5;
parameters.Across_sun = 10;                                % Cross area related to the sun [m^2]
parameters.t0sym = date2mjd2000([2021, 1, 1, 0, 0, 0]);    
parameters.tmax = date2mjd2000([2021, 1, 10, 0, 0, 0]);

parameters.opt = odeset('RelTol',1e-10, 'AbsTol',1e-10,'Events', @event_cont_thrust);

%capture delta_v
kep_cap = [8500 0 0 0 0 0];
delta = deg2rad(175);
[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);
mu = astroConstants(14);

[dt_hyp, dv_req, theta_inf, kep_hyp] = capture(kep_cap, delta,rM, mu);
%the idea is to start the burn at (t_pericenter - t_BO/2) to emulate at
%most an impulsive maneuver at pericenter

