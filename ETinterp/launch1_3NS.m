%% first launch
[kepEarth_NS1, muS]   = uplanet(t01,3);
[kepMars_NS1, ~]    = uplanet(t01 + TOF1,4);

[RE1, vE1] = kep2car2(kepEarth_NS1, muS); %km....
[R2_1, varr_1]   = kep2car2(kepMars_NS1, muS);  %km....

%definition of plane of motion
r1norm_1 = norm(RE1);
r2norm_1 = norm(R2_1);

r1vers_1 = RE1/r1norm_1;
r2vers_1 = R2_1/r2norm_1;

hh_1 = cross(r1vers_1, r2vers_1);
hvers_1 = hh_1/norm(hh_1);

if hh_1(3) <0
    hvers_1 = -hvers_1;
end

%adding TMI maneuver
vE1 = vE1 + v_inf1*(sin(beta1)*cos(alpha1)*r1vers_1 + ...
                 sin(beta1)*sin(alpha1)*cross(hvers_1,r1vers_1) + ...
                 cos(beta1)*hvers_1);

%vinf @ Mars = 0;
             
[ m1, T1, r1, z1, ~, vr1, vt1, vz1, acc_inplane1, acc_out1, acc1, TH1, ~, gamma11, gamma21, gamma1, v1perp1, v2perp1, v1tra1, v2tra1, vnorm1, ~, T_inplane1, T_outplane1, theta_dot1, time1, TOFr1] = ... 
    Conway(TOF1, N_rev1, q1, r1norm_1, r2norm_1, r1vers_1, r2vers_1, hvers_1, hh_1, vE1, varr_1, muS, dataNS);


%% (RS equatorial)

deltaRS1 = 3; %days

TOFRS1 = TOF1 + deltaRS1; 

[kepMars_NS1, ~]    = uplanet(t01 + TOFRS1,4);
[R2_RS1, v2_RS1] = kep2car2(kepMars_NS1, muS);  %km....

%definition of plane of motion
r2norm_RS1 = norm(R2_RS1);
r2vers_RS1 = R2_RS1/r2norm_RS1;

hh_RS1 = cross(r1vers_1, r2vers_RS1);
hvers_RS1 = hh_RS1/norm(hh_RS1);

if hh_RS1(3) <0
    hvers_RS1 = -hvers_RS1;
end

%vinf @ Mars = 0;
             
[ mRS1, TRS1, rRS1, zRS1, ~, vrRS1, vtRS1, vzRS1, acc_inplaneRS1, acc_outRS1, acc_totRS1, THRS1, ~, gamma1RS1, gamma2RS1, gammaRS1, v1perpRS1, v2perpRS1, v1traRS1, v2traRS1, vnormRS1, ~, T_inplaneRS1, T_outplaneRS1, theta_dot1, timeRS1, TOFrRS1] = ... 
    Conway(TOFRS1, N_rev1, q1, r1norm_1, r2norm_RS1, r1vers_1, r2vers_RS1, hvers_RS1, hh_RS1, vE1, v2_RS1, muS, dataRS1);

%% output
fprintf('\n %%%%%% FIRST LAUNCH (NS)\n \n')
[rrNS1, vvNS1, TNS1] = conway3Dplots(t01, TOF1, hvers_1, N_rev1, q1, m1, T1, r1, z1, vr1, vt1, vz1, acc_inplane1, acc_out1, acc1, TH1, gamma11, gamma21, gamma1, v1perp1, v2perp1, v1tra1, v2tra1, vnorm1, T_inplane1, T_outplane1, TOFr1, dataNS);

fprintf('\n %%%%%%  (RS equatorial) \n \n')
[rrRS1, vvRS1, TRS1cart] = conway3Dplots(t01, TOFRS1, hvers_RS1, N_rev1, q1, mRS1, TRS1, rRS1, zRS1, vrRS1, vtRS1, vzRS1, acc_inplaneRS1, acc_outRS1, acc_totRS1, THRS1, gamma1RS1, gamma2RS1, gammaRS1, v1perpRS1, v2perpRS1, v1traRS1, v2traRS1, vnormRS1, T_inplaneRS1, T_outplaneRS1, TOFrRS1, dataRS1);

%% struct
dataNS1.t0MJD2000           = t01;
dataNS1.comps.gamma         = gamma1;
dataNS1.comps.Tacc_in       = acc_inplane1;
dataNS1.comps.Tacc_out      = acc_out1;
dataNS1.comps.T_in          = T_inplane1;
dataNS1.comps.T_out         = T_outplane1;
dataNS1.TOFdays             = TOFr1;
dataNS1.m                   = m1;
dataNS1.T_cart              = TNS1;
dataNS1.r_cart              = rrNS1;
dataNS1.v_cart              = vvNS1;

dataRS1.t0MJD2000           = t01;
dataRS1.comps.gamma         = gammaRS1;
dataRS1.comps.Tacc_in       = acc_inplaneRS1;
dataRS1.comps.Tacc_out      = acc_outRS1;
dataRS1.comps.T_in          = T_inplaneRS1;
dataRS1.comps.T_out         = T_outplaneRS1;
dataRS1.TOFdays             = TOFrRS1;
dataRS1.m                   = mRS1;
dataRS1.T_cart              = TRS1cart;
dataRS1.r_cart              = rrRS1;
dataRS1.v_cart              = vvRS1;