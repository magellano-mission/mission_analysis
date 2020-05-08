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
% parameters.c_r = 0.5;
% parameters.Across_sun = 10;                                % Cross area related to the sun [m^2]
  
parameters.event = 0;
parameters.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12);   

mu = astroConstants(14);

%capture delta_v
% kep_cap_desired = [10000 0.1 deg2rad(90) 0 0 0];


%% first launch (2024)
parameters.t0sym = 9.361401106238969e+03;  %MJD2000

[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);

VF =    [12.799958274752937 -16.478235594522850   0.335452507074074]; %SM heliocentric velocity
v_M =   [13.837765883347817 -18.403996797073731  -0.725009406810734]; %mars heliocentric velocity

kep_NS = [9850 0 deg2rad(25) 0 0 0];
kep_NS2 = [8100 0 deg2rad(25) 0 0 0];

[YY_capture, hyp_capture , kep_capture] = PO2hyp(kep_NS, (VF - v_M), rM, mu, parameters, Thrust0, 1, 'arrival', 2);
%%
[YY_capture, hyp_capture , kep_capture] = PO2hyp(kep_NS2, (VF - v_M), rM, mu, parameters, Thrust0, 1, 'arrival', 2);

%% second launch (2026)
parameters.t0sym = 1.010485462521553e+04;  %MJD200

[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);

VF =    [19.599339413692434  -8.951498304070178  -0.215260206190249]; %SM heliocentric velocity
v_M =   [21.963857541381959  -9.885170051123481  -0.747631467874832]; %mars heliocentric velocity

kep_NS3 = [10500 0 deg2rad(25) 0 0 0]; 
kep_RS = [6400 0 deg2rad(1e-8) 0 0 0];

[YY_capture, hyp_capture , kep_capture] = PO2hyp(kep_NS3, (VF - v_M), rM, mu, parameters, Thrust0, 1, 'arrival', 2);
%%
[YY_capture, hyp_capture , kep_capture] = PO2hyp(kep_RS, (VF - v_M), rM, mu, parameters, Thrust0, 1, 'arrival', 2);


