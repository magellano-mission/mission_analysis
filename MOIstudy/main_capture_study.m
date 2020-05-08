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
parameters.isPerturbed = 0;                                % 

% parameters.dt_p = parameters.t_BO/2;  %defined inside                      %time distance from pericenter at which firing occurs [s]
% parameters.t_BO = 30*60;   %computed inside              % burnout time (chemical thrust)
parameters.T = [0; 0; 0];                                  % thrust [N] (@TNH) (constant profile for the moment)
parameters.perc_tBO_before_p = 0.5;
Thrust0 = [-10000; 0; 0];                                   % thrust [N] (@TNH) (constant profile for the moment)
parameters.Isp = 380;                                      % specific impulse [s]

parameters.M0 = 6000;                                       % Total Mass of the s/c [kg]
parameters.c_r = 0.5;
parameters.Across_sun = 10;                                % Cross area related to the sun [m^2]
parameters.t0sym = date2mjd2000([2021, 1, 1, 0, 0, 0]);    
parameters.event = 0;
parameters.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12);   %,'Events', @event_cont_thrust);


%capture delta_v
kep_cap_desired = [10000 0.1 deg2rad(7) 0 0 0];
delta = deg2rad(45);
[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);
mu = astroConstants(14);

VF =    [20.1823   -7.8639   -0.3019]; %SM heliocentric velocity
v_M =   [22.5170   -8.8545   -0.7398]; %mars heliocentric velocity

[YY_capture, hyp_capture , kep_capture] = PO2hyp(kep_cap_desired, (VF - v_M), rM, mu, parameters, Thrust0, 1, 'arrival', 2);


