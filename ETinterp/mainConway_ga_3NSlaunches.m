%Conway-based approach
run('MatlabSetup.m')
%% SPACECRAFT CHARACTERISTICS
data.Isp = 4300;                                       % specific impulse [s]
data.Mdry = 1600;                                      % Total Mass of the s/c [kg]
data.n_int = 1000;
data.Tmax = 0.2; % N
data.n_int = 1000;

choice = menu('Do you want to run a new simulation or load a SOL?','Run a simulation','Load SOL');
if choice == 1 || choice == 0
%% BOUNDARIES
% t0, TOF, N_rev, q, v_inf, alpha, beta, v_infcap, alphacap, betacap 
%%%%%
% lower boundary
lb = zeros(1,21); ub = lb;
lb(1) = date2mjd2000([2024 1 1 0 0 0]); %t0
lb(2) = 900; %TOF
lb(3) = 1; %N_rev
lb(4) = 3;% q
lb(5) = 0;% v_inf dep
lb(6) = 0;% alpha dep
lb(7) = 0;% beta dep
lb(8) = 900; %TOF
lb(9) = 1; %N_rev
lb(10) = 3;% q
lb(11) = 0;% v_inf dep
lb(12) = 0;% alpha dep
lb(13) = 0;% beta dep
lb(14) = 1000; %TOF
lb(15) = 1; %N_rev
lb(16) = 3;% q
lb(17) = 0;% v_inf dep
lb(18) = 0;% alpha dep
lb(19) = 0;% beta dep

% upper boundary 
ub(1) = date2mjd2000([2026 1 1 0 0 0]);
ub(2) = 1000;
ub(3) = 3;
ub(4) = 8;
ub(5) = 0.5; %v_inf dep
ub(6) = pi;% alpha dep
ub(7) = pi;% beta dep
ub(8) = 1000;
ub(9) = 3;
ub(10) = 8;
ub(11) = 1.5; %v_inf dep
ub(12) = pi;% alpha dep
ub(13) = pi;% beta dep
ub(14) = 1100;
ub(15) = 3;
ub(16) = 8;
ub(17) = 1.5; %v_inf dep
ub(18) = pi;% alpha dep
ub(19) = pi;% beta dep

%flexibility margins
lb(20) = -5;
lb(21) = -5;

ub(20) = 5;
ub(21) = 5;
 
Bound = [lb; ub];

%% choose optimization method
optimization_method = menu('Ok! Remember to set boudaries... \n Choose the optimization method','ga','gamultiboj','particleswarm');


switch optimization_method
    case 1
    options = optimoptions('ga', 'Display', 'Iter', ...
                           'PopulationSize', 200, 'StallGenLimit', 200, ...         
                           'MaxGenerations', 200, ...
                           'UseParallel', true, 'PopInitRange',Bound);
    [SOL,feval,exitflag] = ga(@(x) ga_conway_3NSlaunches(x,data), 21,[],[],[],[],lb,ub,[], [3 9 15],options);
    chosen = 1;
    case 2
    options = optimoptions('gamultiobj', 'Display', 'Iter', ...
                           'PopulationSize', 200, 'StallGenLimit', 200, ...         
                           'MaxGenerations', 200, ...
                           'ParetoFraction', 0.35, ...
                           'UseParallel', true, 'PopInitRange',Bound);
    [SOL,feval,exitflag] = gamultiobj(@(x) gamultiobj_conway_3NSlaunches(x,data), 21,[],[],[],[],lb,ub,[],options);
    [chosen, chosenT] = paretoplot(SOL, feval,data); %min TOF chosen
    chosen = chosenT;
    case 3
    options = optimoptions('particleswarm', 'SwarmSize', 1000, 'MaxIterations', 20, 'Display', 'Iter','UseParallel', true);
    [SOL,feval,exitflag] = particleswarm(@(x) ga_conway_3NSlaunches(x, data), 21, lb, ub, options);
    chosen = 1;
    case 0
        fprintf('\n loaded a solution... \n')
    load('SOLbest.mat')
    chosen = 1;
end
else
    load('SOLbest.mat')
    chosen = 1;
end
clear choice

%% SOLUTIONS
t01 = SOL(chosen, 1);             
TOF1 = SOL(chosen, 2);
N_rev1 = round(SOL(chosen, 3));    
q1 = SOL(chosen, 4);
v_inf1 = SOL(chosen, 5);
alpha1 = SOL(chosen, 6);
beta1 = SOL(chosen, 7);

%second launch
TOF2 = SOL(chosen, 8);
N_rev2 = round(SOL(chosen, 9));  
q2 = SOL(chosen, 10);
v_inf2 = SOL(chosen, 11);
alpha2 = SOL(chosen, 12); 
beta2 = SOL(chosen, 13);

%flexiblity 
delta1 = SOL(chosen,20);
t02 = t01 + TOF1 + 186 - TOF2 + delta1;

%third launch
TOF3 = SOL(chosen, 14);
N_rev3 = round(SOL(chosen, 15));    
q3 = SOL(chosen, 16);
v_inf3 = SOL(chosen, 17); 
alpha3 = SOL(chosen, 18); 
beta3 = SOL(chosen, 19);

delta2 = SOL(chosen,21);
t03 = t01 + TOF1 + 372 - TOF3 + delta2;

%first launch
[kepEarth_1, muS] = uplanet(t01,3);
[kepMars_1, ~]    = uplanet(t01 + TOF1,4);

[R1_1, v1_1] = kep2car2(kepEarth_1, muS); %km....
[R2_1, v2_1] = kep2car2(kepMars_1, muS);  %km....

%definition of plane of motion
r1norm_1 = norm(R1_1);
r2norm_1 = norm(R2_1);

r1vers_1 = R1_1/r1norm_1;
r2vers_1 = R2_1/r2norm_1;

hh_1 = cross(r1vers_1, r2vers_1);
hvers_1 = hh_1/norm(hh_1);

if hh_1(3) <0
    hvers_1 = -hvers_1;
end

%adding TMI maneuver
v1_1 = v1_1 + v_inf1*(sin(beta1)*cos(alpha1)*r1vers_1 + ...
                 sin(beta1)*sin(alpha1)*cross(hvers_1,r1vers_1) + ...
                 cos(beta1)*hvers_1);

%vinf @ Mars = 0;
             
[ m1, T1, r1, z1, ~, vr1, vt1, vz1, acc_inplane1, acc_out1, acc1, TH1, ~, gamma11, gamma21, gamma1, v1perp1, v2perp1, v1tra1, v2tra1, vnorm1, ~, T_inplane1, T_outplane1, theta_dot1, time1, TOFr1] = ... 
    Conway(TOF1, N_rev1, q1, r1norm_1, r2norm_1, r1vers_1, r2vers_1, hvers_1, hh_1, v1_1, v2_1, muS, data);


%second launch
[kepEarth_2, ~] = uplanet(t02,3);
[kepMars_2, ~]    = uplanet(t02 + TOF2,4);

[R1_2, v1_2] = kep2car2(kepEarth_2, muS); %km....
[R2_2, v2_2] = kep2car2(kepMars_2, muS);  %km....

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
    Conway(TOF2, N_rev2, q2, r1norm_2, r2norm_2, r1vers_2, r2vers_2, hvers_2, hh_2, v1_2, v2_2, muS, data);

%thrid launch
[kepEarth_3, ~] = uplanet(t03,3);
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
    Conway(TOF3, N_rev3, q3, r1norm_3, r2norm_3, r1vers_3, r2vers_3, hvers_3, hh_3, v1_3, v2_3, muS, data);

fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%%%')
fprintf('\n')
fprintf('\n %%%%%% FIRST LAUNCH \n \n')
[rrNS1, vvNS1] = conway3Dplots(t01, TOF1, hvers_1, N_rev1, q1, m1, T1, r1, z1, vr1, vt1, vz1, acc_inplane1, acc_out1, acc1, TH1, gamma11, gamma21, gamma1, v1perp1, v2perp1, v1tra1, v2tra1, vnorm1, T_inplane1, T_outplane1, TOFr1, data);

fprintf('\n press ENTER to continue... \n')
pause()

fprintf('\n %%%%%% SECOND LAUNCH \n \n')
[rrNS2, vvNS2] = conway3Dplots(t02, TOF2, hvers_2, N_rev2, q2, m2, T2, r2, z2, vr2, vt2, vz2, acc_inplane2, acc_out2, acc2, TH2, gamma12, gamma22, gamma2, v1perp2, v2perp2, v1tra2, v2tra2, vnorm2, T_inplane2, T_outplane2, TOFr2, data);

fprintf('\n press ENTER to continue... \n')
pause()

fprintf('\n %%%%%% THIRD LAUNCH \n \n')
[rrNS3, vvNS3] = conway3Dplots(t03, TOF3, hvers_3, N_rev3, q3, m3, T3, r3, z3, vr3, vt3, vz3, acc_inplane3, acc_out3, acc3, TH3, gamma13, gamma23, gamma3, v1perp3, v2perp3, v1tra3, v2tra3, vnorm3, T_inplane3, T_outplane3, TOFr3, data);

fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%%%')
fprintf('\n')

%%
% MONTECARLO SIMULATION
%giving the same control law (no predictive control) what happens to the
%state...
NMC = 100;
X0 = [r1(1), TH1(1), z1(1), vr1(1), theta_dot1(1), vz1(1), m1(1)];
XMC = zeros(data.n_int, 7,NMC); 
sigma_inj = X0*0.001.*ones(1,7); 
sigma_inj(7) = 0; %same mass of propellant of course...

wb = waitbar(0,'MONTECARLO SIMULATION');
for nmc = 1:NMC
   waitbar(nmc/NMC, wb);
   X0MC = X0 + sigma_inj.*randn(1,7);
  [XMC(:,:,nmc)] = EoMpropRK4(X0MC, time1, acc_inplane1, acc_out1, muS, data, 0);
end   
delete wb
rrMC = zeros(data.n_int, 3, NMC);vvMC = rrMC;
pause()
%computation of the state of 
for nmc = 1:NMC
    r   = squeeze(XMC(:,1,nmc));
    th  = squeeze(XMC(:,2,nmc));
    z   = squeeze(XMC(:,3,nmc));
    vr  = squeeze(XMC(:,4,nmc));
    vt  = r.*squeeze(XMC(:,5,nmc));
    vz  = squeeze(XMC(:,6,nmc));
    for i = 1:data.n_int
        [rrMC(i,:,nmc), vvMC(i,:,nmc)] = refplane2car( r(i), z(i),  vt(i), vr(i), vz(i), th(i), r1vers_1, hvers_1);
    end
end
%%
figure()
RE = zeros(length(TOFr1), 3); RM = RE;
for i =1:length(TOFr1)
    kepE = uplanet(t01+TOFr1(i), 3);
    kepM = uplanet(t01+TOFr1(i), 4);
    RE(i,:) = kep2car2(kepE, muS);
    RM(i,:) = kep2car2(kepM, muS);
end
clear kepE kepM
mpp = plot3(RM(:,1),RM(:,2),RM(:,3),'HandleVisibility','Off'); hold on
epp = plot3(RE(:,1),RE(:,2),RE(:,3),'HandleVisibility','Off'); hold on
plot3(rrNS1(:,1), rrNS1(:,2), rrNS1(:,3),'m', 'LineWidth',5, 'DisplayName','Conway Solution'), hold on
sgtitle('Injection maneuver error')
for nmc = 1:NMC
    ppp = plot3(rrMC(:,1,nmc),rrMC(:,2,nmc),rrMC(:,3,nmc),'k:','HandleVisibility','Off'); hold on
    ppp.Color(4) = 0.2;
end
ppp.HandleVisibility = 'on';
ppp.DisplayName = 'Montecarlo Simulation';

mpp0 = plot3(RM(end,1),RM(end,2),RM(end,3), 'o','MarkerSize',10,'DisplayName','Mars @ Arrival'); hold on
epp0 = plot3(RE(1,1),RE(1,2),RE(1,3),'o','MarkerSize',10,'DisplayName','Earth @ Departure'); hold on
mpp0.Color = mpp.Color;
epp0.Color = epp.Color;
plotSun(), legend()

function [chosen, chosenT] = paretoplot(SOL, feval, data)
figure()
sgtitle('Pareto Front')
minTOF = 1e5; chosen = 0;
for i = 1:length(feval)
    if feval(i,1) <= data.Tmax
        plot(feval(i,1), feval(i,2), 'ro','HandleVisibility','off'), hold on
    else 
        plot(feval(i,1), feval(i,2), 'k+','HandleVisibility','off'), hold on
    end
         if SOL(i,2) + 372 < minTOF
            minTOF = SOL(i,2) + 372;
            chosen = i;
         end
end
minmaxT1 = feval(find(feval(:,1) == min(feval(:,1))),1) ;
minmaxT2 = feval(find(feval(:,2) == min(feval(:,2))),2) ;
minmaxT3 = feval(find(feval(:,3) == min(feval(:,3))),3) ;
if minmaxT1 > minmaxT2 && minmaxT1 > minmaxT3
    chosenT = find(feval(:,1) == min(feval(:,1)));
end
if minmaxT2 > minmaxT1 && minmaxT2 > minmaxT3
    chosenT = find(feval(:,2) == min(feval(:,2)));
end
if minmaxT3 > minmaxT1 && minmaxT3 > minmaxT2
    chosenT = find(feval(:,3) == min(feval(:,3)));
end

plot(feval(chosen,1), feval(chosen,2),'go', 'DisplayName',strcat('Min set-up time (', num2str(SOL(chosen,2)),'days )')); hold on
plot(feval(chosenT,1), feval(chosenT,2),'bo', 'DisplayName',strcat('Min ( max(T(t))) (', num2str(feval(chosenT,1)),'N )')); hold off
xline(data.Tmax,'r','DisplayName','Max Thrust < 0.2 N'), xlabel('max(|T(t)|) [N]'), ylabel('m_P [kg]') , legend()
end
