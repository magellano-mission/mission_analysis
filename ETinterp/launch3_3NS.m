%% thrid launch
[kepEarth_3, ~]   = uplanet(t03,3);
[kepMars_3, ~]    = uplanet(t03 + TOF3,4);

[R1_3, v1_3] = kep2car2(kepEarth_3, muS); %km....
[R2_3, v2_3] = kep2car2(kepMars_3, muS);  %km....

%definition of plane of motion
r1norm_3 = norm(R1_3);
r2norm_3 = norm(R2_3);

r1vers_3 = R1_3/r1norm_3;
r2vers_3 = R2_3/r2norm_3;

hh_3 = cross(r1vers_3, r2vers_3);
hvers_3 = hh_3/norm(hh_3);

if hh_3(3) <0
    hvers_3 = -hvers_3;
end

%adding TMI maneuver
v1_3 = v1_3 + v_inf3*(sin(beta3)*cos(alpha3)*r1vers_3 + ...
                 sin(beta3)*sin(alpha3)*cross(hvers_3,r1vers_3) + ...
                 cos(beta3)*hvers_3);

%vinf @ Mars = 0;
             
[ m3, T3, r3, z3, ~, vr3, vt3, vz3, acc_inplane3, acc_out3, acc3, TH3, ~, gamma13, gamma23, gamma3, v1perp3, v2perp3, v1tra3, v2tra3, vnorm3, ~, T_inplane3, T_outplane3, theta_dot3, ~, TOFr3] = ...
    Conway(TOF3, N_rev3, q3, r1norm_3, r2norm_3, r1vers_3, r2vers_3, hvers_3, hh_3, v1_3, v2_3, muS, dataNS);

%% output
fprintf('\n %%%%%% THIRD LAUNCH (NS) \n \n')
[rrNS3, vvNS3, TNS3] = conway3Dplots(t03, TOF3, hvers_3, N_rev3, q3, m3, T3, r3, z3, vr3, vt3, vz3, acc_inplane3, acc_out3, acc3, TH3, gamma13, gamma23, gamma3, v1perp3, v2perp3, v1tra3, v2tra3, vnorm3, T_inplane3, T_outplane3, TOFr3, dataNS);

%% struct
dataNS3.t0MJD2000           = t03;
dataNS3.comps.gamma         = gamma3;
dataNS3.comps.Tacc_in       = acc_inplane3;
dataNS3.comps.Tacc_out      = acc_out3;
dataNS3.comps.T_in          = T_inplane3;
dataNS3.comps.T_out         = T_outplane3;
dataNS3.TOFdays             = TOFr3;
dataNS3.m                   = m3;
dataNS3.T_cart              = TNS3;
dataNS3.r_cart              = rrNS3;
dataNS3.v_cart              = vvNS3;

