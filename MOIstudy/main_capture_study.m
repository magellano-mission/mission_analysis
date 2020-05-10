%% CAPTURE STUDY
%Propagation of a single-burn maneuver for a generic capture (more refined
%model coming in one day)
%change the Thrust vector, M0, Isp, (Mp0 is an output at the moment,
%myabe it is better to understand the order of magnitudes involved)
%delta_t_before_pericenter: trying to understand the impact of the firing
%timing

close all, clear, clc

parameters.isInterp = 0;
parameters.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
parameters.isPerturbed = 0;                                % 

parameters.T = [0; 0; 0];                                  % thrust [N] (@TNH) (constant profile for the moment)
parameters.perc_tBO_before_p = 0.5;
Thrust0 = [-420*5; 0; 0];                                   % thrust [N] (@TNH) (constant profile for the moment)
parameters.Isp = 320;                                      % specific impulse [s]

parameters.M0 = 7000;                                       % Total Mass of the s/c [kg]
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

[YYNS, hyp_capture , kep_capture] = hyp2PO(kep_NS, (VF - v_M), rM, mu, parameters, Thrust0, 1, 'arrival', 2);
[YYNS2, hyp_capture , kep_capture] = hyp2PO(kep_NS2, (VF - v_M), rM, mu, parameters, Thrust0, 1, 'arrival', 2);

figure()
I = imread('Mars.jpg'); RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];  RI.YWorldLimits = [-90 90]; 
rMars = almanac('Mars','Radius','kilometers','sphere');
[X, Y, Z] = ellipsoid(0, 0, 0, rMars, rMars, rMars, 100); % spheric centered Mars
planet = surf(X, Y, -Z,'Edgecolor', 'none','DisplayName', 'Sun'); hold on
set(planet,'FaceColor','texturemap','Cdata',I), axis equal

cap(1) = plot3(YYNS(:,1), YYNS(:,2), YYNS(:,3), 'DisplayName','NS');
cap(2) = plot3(YYNS2(:,1), YYNS2(:,2), YYNS2(:,3), 'DisplayName','NS2');
xlim([-11000 11000]), ylim([-11000 11000]), zlim([-11000 11000])

cap(1).Color = [0.9490    0.4745    0.3137];
cap(2).Color = [0.1020    0.6667    0.74120];

legend()
%% second launch (2026)
parameters.t0sym = 1.010485462521553e+04;  %MJD200

[kepM, muS] = uplanet(parameters.t0sym, 4);
rM = kep2car2(kepM, muS);

VF =    [19.599339413692434  -8.951498304070178  -0.215260206190249]; %SM heliocentric velocity
v_M =   [21.963857541381959  -9.885170051123481  -0.747631467874832]; %mars heliocentric velocity

kep_NS3 = [10500 0 deg2rad(25) 0 0 0]; 
kep_RS = [6400 0 deg2rad(1e-8) 0 0 0];

[YYNS3, hyp_capture , kep_capture] = hyp2PO(kep_NS3, (VF - v_M), rM, mu, parameters, Thrust0, 1, 'arrival', 2);
[YYRS, hyp_capture , kep_capture] = hyp2PO(kep_RS, (VF - v_M), rM, mu, parameters, Thrust0, 1, 'arrival', 2);

figure()
I = imread('Mars.jpg'); RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];  RI.YWorldLimits = [-90 90]; 
rMars = almanac('Mars','Radius','kilometers','sphere');
[X, Y, Z] = ellipsoid(0, 0, 0, rMars, rMars, rMars, 100); % spheric centered Mars
planet = surf(X, Y, -Z,'Edgecolor', 'none','DisplayName', 'Sun'); hold on
set(planet,'FaceColor','texturemap','Cdata',I), axis equal

plot3(YYNS(:,1), YYNS(:,2), YYNS(:,3), 'k--', 'DisplayName','NS');
plot3(YYNS2(:,1), YYNS2(:,2), YYNS2(:,3), 'k--', 'DisplayName','NS2');

cap(1) = plot3(YYNS3(:,1), YYNS3(:,2), YYNS3(:,3), 'DisplayName','NS + ECS');
cap(2) = plot3(YYRS(:,1), YYRS(:,2), YYRS(:,3), 'DisplayName','RS');
xlim([-11000 11000]), ylim([-11000 11000]), zlim([-11000 11000])

cap(1).Color = [0.9490    0.4745    0.3137];
cap(2).Color = [0.1020    0.6667    0.74120];

legend()
