%% second launch
[kepEarth_2, ~] = uplanet(t02,3);
[kepMars_RS2, ~]    = uplanet(t02 + TOF2,4);

[R1_2, v1_2] = kep2car2(kepEarth_2, muS); %km....
[R2_2, v2_2] = kep2car2(kepMars_RS2, muS);  %km....

%definition of plane of motion
r1norm_2 = norm(R1_2);
r2norm_2 = norm(R2_2);

r1vers_2 = R1_2/r1norm_2;
r2vers_2 = R2_2/r2norm_2;

hh_2 = cross(r1vers_2, r2vers_2);
hvers_2 = hh_2/norm(hh_2);

if hh_2(3) <0
    hvers_2 = -hvers_2;
end

%adding TMI maneuver
v1_2 = v1_2 + v_inf2*(sin(beta2)*cos(alpha2)*r1vers_2 + ...
                 sin(beta2)*sin(alpha2)*cross(hvers_2,r1vers_2) + ...
                 cos(beta2)*hvers_2);

%vinf @ Mars = 0;
             
[ m2, T2, r2, z2, ~, vr2, vt2, vz2, acc_inplane2, acc_out2, acc2, TH2, ~, gamma12, gamma22, gamma2, v1perp2, v2perp2, v1tra2, v2tra2, vnorm2, ~, T_inplane2, T_outplane2, theta_dot2, ~, TOFr2] = ...
    Conway(TOF2, N_rev2, q2, r1norm_2, r2norm_2, r1vers_2, r2vers_2, hvers_2, hh_2, v1_2, v2_2, muS, dataNS);

%% (RS polar)

deltaRS2 = 3; %days

TOFRS2 = TOF2 + deltaRS2; 

[kepMars_RS2, ~]    = uplanet(t02 + TOFRS2,4);
[R2_RS2, v2_RS2] = kep2car2(kepMars_RS2, muS);  %km....

%definition of plane of motion
r2norm_RS2 = norm(R2_RS2);
r2vers_RS2 = R2_RS2/r2norm_RS2;

hh_RS2 = cross(r1vers_2, r2vers_RS2);
hvers_RS2 = hh_RS2/norm(hh_RS2);

if hh_RS2(3) <0
    hvers_RS2 = -hvers_RS2;
end

%vinf @ Mars = 0;
             
[ mRS2, TRS2, rRS2, zRS2, ~, vrRS2, vtRS2, vzRS2, acc_inplaneRS2, acc_outRS2, acc_totRS2, THRS2, ~, gamma1RS2, gamma2RS2, gammaRS2, v1perpRS2, v2perpRS2, v1traRS2, v2traRS2, vnormRS2, ~, T_inplaneRS2, T_outplaneRS2, theta_dot2, timeRS2, TOFrRS2] = ... 
    Conway(TOFRS2, N_rev2, q2, r1norm_2, r2norm_RS2, r1vers_2, r2vers_RS2, hvers_RS2, hh_RS2, v1_2, v2_RS2, muS, dataRS2);

%% output
fprintf('\n %%%%%% SECOND LAUNCH (NS) \n \n')
[rrNS2, vvNS2, TNS2] = conway3Dplots(t02, TOF2, hvers_2, N_rev2, q2, m2, T2, r2, z2, vr2, vt2, vz2, acc_inplane2, acc_out2, acc2, TH2, gamma12, gamma22, gamma2, v1perp2, v2perp2, v1tra2, v2tra2, vnorm2, T_inplane2, T_outplane2, TOFr2, dataNS);

fprintf('\n %%%%%%  (RS polar) \n \n')
[rrRS2, vvRS2, TRS2cart] = conway3Dplots(t02, TOFRS2, hvers_RS2, N_rev2, q2, mRS2, TRS2, rRS2, zRS2, vrRS2, vtRS2, vzRS2, acc_inplaneRS2, acc_outRS2, acc_totRS2, THRS2, gamma1RS2, gamma2RS2, gammaRS2, v1perpRS2, v2perpRS2, v1traRS2, v2traRS2, vnormRS2, T_inplaneRS2, T_outplaneRS2, TOFrRS2, dataRS2);

%% struct
dataNS2.t0MJD2000           = t02;
dataNS2.comps.gamma         = gamma2;
dataNS2.comps.Tacc_in       = acc_inplane2;
dataNS2.comps.Tacc_out      = acc_out2;
dataNS2.comps.T_in          = T_inplane2;
dataNS2.comps.T_out         = T_outplane2;
dataNS2.TOFdays             = TOFr2;
dataNS2.m                   = m2;
dataNS2.T_cart              = TNS2;
dataNS2.r_cart              = rrNS2;
dataNS2.v_cart              = vvNS2;

dataRS2.t0MJD2000           = t02;
dataRS2.comps.gamma         = gammaRS2;
dataRS2.comps.Tacc_in       = acc_inplaneRS2;
dataRS2.comps.Tacc_out      = acc_outRS2;
dataRS2.comps.T_in          = T_inplaneRS2;
dataRS2.comps.T_out         = T_outplaneRS2;
dataRS2.TOFdays             = TOFrRS2;
dataRS2.m                   = mRS2;
dataRS2.T_cart              = TRS2cart;
dataRS2.r_cart              = rrRS2;
dataRS2.v_cart              = vvRS2;