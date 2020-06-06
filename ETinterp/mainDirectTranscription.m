%% Direct Transcription
% FIRST LAUNCH
dataNS1.n_int = 70; 

Bounds.alpha_lb = deg2rad(-10);%-pi;
Bounds.alpha_ub = deg2rad(10);%pi;
Bounds.beta_lb = deg2rad(-20);%-pi/4;
Bounds.beta_ub = deg2rad(20);%pi/4;

%r
Bounds.r_ub = 1.55;
Bounds.r_lb =  0.85;

%theta
Bounds.th_ub = Inf;
Bounds.th_lb = 0;
%z
Bounds.z_ub =   0.05;
Bounds.z_lb =  -0.05;
%vr
Bounds.vr_ub = 0.1;
Bounds.vr_lb = -0.2;
%theta_dot
Bounds.thd_ub = Inf;
Bounds.thd_lb = -Inf;
%vz
Bounds.vz_ub = 0.03;
Bounds.vz_lb = -0.03;

[THSNS1, alphaHSNS1, betaHSNS1, XHSNS1, XSOLNS1,fvalNS1,exitflagNS1,outputNS1,lambdaNS1,gradNS1,hessianNS1] = ...
runDirectTranscription(t01, TOF1, N_rev1, q1, r1norm_1, r2norm_1, r1vers_1, r2vers_1, hvers_1, hh_1, v1_1, v2_1, muS, dataNS1, Bounds);
%%
dataRS1.n_int = 70;

Bounds.alpha_lb = deg2rad(-10);%-pi;
Bounds.alpha_ub = deg2rad(10);%pi;
Bounds.beta_lb = deg2rad(-45);%-pi/4;
Bounds.beta_ub = deg2rad(45);%pi/4;

%r
Bounds.r_ub = 1.55;
Bounds.r_lb = 0.85;

%theta
Bounds.th_ub = Inf;
Bounds.th_lb = 0;
%z
Bounds.z_ub =  0.05;
Bounds.z_lb =  -0.05;
%vr
Bounds.vr_ub = 0.15;
Bounds.vr_lb = -0.15;
%theta_dot
Bounds.thd_ub = Inf;
Bounds.thd_lb = -Inf;
%vz
Bounds.vz_ub = 0.05;
Bounds.vz_lb = -0.05;
    
[THSRS1, alphaHSRS1, betaHSRS1, XHSRS1, XSOLRS1,fvalRS1,exitflagRS1,outputRS1,lambdaRS1,gradRS1,hessianRS1] = ...
runDirectTranscription(t01, TOFRS1, N_rev1, q1, r1norm_1, r2norm_RS1, r1vers_1, r2vers_RS1, hvers_RS1, hh_RS1, v1_1, v2_RS1, muS, dataRS1, Bounds);

%%
dataECS.n_int = 70; 

Bounds.alpha_lb = deg2rad(-45);%-pi;
Bounds.alpha_ub = deg2rad(45);%pi;
Bounds.beta_lb = deg2rad(-45);%-pi/4;
Bounds.beta_ub = deg2rad(45);%pi/4;

%r
Bounds.r_ub = 1.55;
Bounds.r_lb = 0.85;

%theta
Bounds.th_ub = Inf;
Bounds.th_lb = 0;
%z
Bounds.z_ub =  0.05;
Bounds.z_lb =  -0.05;
%vr
Bounds.vr_ub = 0.15;
Bounds.vr_lb = -0.15;
%theta_dot
Bounds.thd_ub = Inf;
Bounds.thd_lb = -Inf;
%vz
Bounds.vz_ub = 0.05;
Bounds.vz_lb = -0.05;

[THSECS, alphaHSECS, betaHSECS, XHSECS, XSOLECS,fvalECS,exitflagECS,outputECS,lambdaECS,gradECS,hessianECS] = ...
runDirectTranscription(t01, TOFECS, N_rev1, q1, r1norm_1, r2norm_ECS, r1vers_1, r2vers_ECS, hvers_ECS, hh_ECS, v1_1, v2_ECS, muS, dataECS, Bounds);
%% SECOND LAUNCH
dataNS2.n_int = 70;

Bounds.alpha_lb = deg2rad(-45);%-pi;
Bounds.alpha_ub = deg2rad(45);%pi;
Bounds.beta_lb = deg2rad(-45);%-pi/4;
Bounds.beta_ub = deg2rad(45);%pi/4;
%r
Bounds.r_ub = 1.80;
Bounds.r_lb = 0.85;

%theta
Bounds.th_ub = Inf;
Bounds.th_lb = 0;
%z
Bounds.z_ub =  0.05;
Bounds.z_lb =  -0.05;
%vr
Bounds.vr_ub = 0.15;
Bounds.vr_lb = -0.15;
%theta_dot
Bounds.thd_ub = Inf;
Bounds.thd_lb = -Inf;
%vz
Bounds.vz_ub = 0.05;
Bounds.vz_lb = -0.05;

[THSNS2, alphaHSNS2, betaHSNS2, XHSNS2, XSOLNS2,fvalNS2,exitflagNS2,outputNS2,lambdaNS2,gradNS2,hessianNS2] = ...
runDirectTranscription (t02, TOF2, N_rev2, q2, r1norm_2, r2norm_2, r1vers_2, r2vers_2, hvers_2, hh_2, v1_2, v2_2, muS, dataNS2, Bounds);

%%
dataRS2.n_int = 70;

Bounds.alpha_lb = deg2rad(-45);%-pi;
Bounds.alpha_ub = deg2rad(45);%pi;
Bounds.beta_lb = deg2rad(-45);%-pi/4;
Bounds.beta_ub = deg2rad(45);%pi/4;
%r
Bounds.r_ub = 1.80;
Bounds.r_lb = 0.85;

%theta
Bounds.th_ub = Inf;
Bounds.th_lb = 0;
%z
Bounds.z_ub =  0.05;
Bounds.z_lb =  -0.05;
%vr
Bounds.vr_ub = 0.15;
Bounds.vr_lb = -0.15;
%theta_dot
Bounds.thd_ub = Inf;
Bounds.thd_lb = -Inf;
%vz
Bounds.vz_ub = 0.05;
Bounds.vz_lb = -0.05;

[THSRS2, alphaHSRS2, betaHSRS2, XHSRS2, XSOLRS2,fvalRS2,exitflagRS2,outputRS2,lambdaRS2,gradRS2,hessianRS2] = ...
runDirectTranscription(t02, TOFRS2, N_rev2, q2, r1norm_2, r2norm_RS2, r1vers_2, r2vers_RS2, hvers_RS2, hh_RS2, v1_2, v2_RS2, muS, dataRS2, Bounds);
%%
dataNS3.n_int = 70;

Bounds.alpha_lb = deg2rad(-35);%-pi;
Bounds.alpha_ub = deg2rad(35);%pi;
Bounds.beta_lb = deg2rad(-45);%-pi/4;
Bounds.beta_ub = deg2rad(45);%pi/4;
%r
Bounds.r_ub = 1.85;
Bounds.r_lb = 0.85;

%theta
Bounds.th_ub = Inf;
Bounds.th_lb = 0;
%z
Bounds.z_ub =  0.05;
Bounds.z_lb =  -0.05;
%vr
Bounds.vr_ub = 0.15;
Bounds.vr_lb = -0.15;
%theta_dot
Bounds.thd_ub = Inf;
Bounds.thd_lb = -Inf;
%vz
Bounds.vz_ub = 0.05;
Bounds.vz_lb = -0.05;

[THSNS3, alphaHSNS3, betaHSNS3, XHSNS3, XSOLNS3,fvalNS3,exitflagNS3,outputNS3,lambdaNS3,gradNS3,hessianNS3] = ...
runDirectTranscription(t03, TOF3, N_rev3, q3, r1norm_3, r2norm_3, r1vers_3, r2vers_3, hvers_3, hh_3, v1_3, v2_3, muS, dataNS3, Bounds);
