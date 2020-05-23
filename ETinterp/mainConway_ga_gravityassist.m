%Conway-based with one gravity assist

%definition of parameters
close all, clear, clc
% Figure Initialization    
load('MagellanoColorMap.mat');
DefaultOrderColor = get(0, 'DefaultAxesColorOrder');
NewOrderColor = [0.9490    0.4745    0.3137
                 0.1020    0.6667    0.74120
                 155/255   155/255   155/255
                 DefaultOrderColor];  
             
set(0,'DefaultFigureColormap', MagellanoColorMap);
set(0, 'DefaultAxesColorOrder', NewOrderColor);
set(0,'DefaultLineLineWidth', 2)
set(0,'DefaultLineMarkerSize', 10)
set(0, 'DefaultFigureUnits', 'normalized');
set(0, 'DefaultFigurePosition', [0 0 1 1]);
set(0, 'DefaultTextFontSize', 18);
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

data_stacks.Isp = 3000;                                      % specific impulse [s]
data_stacks.Mdry = 7000;                                      % Total Mass of the s/c [kg]
data_stacks.n_int = 1000;

%%%%%
% t0 % TOF12 % TOF23 
% N_rev1 % N_rev2 % q1 % q2 
% v_inf % alpha % beta 
% v_infga % alphagam % alphagap % betagam % betagap 
%%%%%
%lower boundary
lb = zeros(1,15); ub = lb;
lb(1) = date2mjd2000([2024 1 1 0 0 0]);  % t0    
lb(2) = 200; % TOF12   
lb(3) = 700; % TOF23   
lb(4) = 0; % N_rev1
lb(5) = 0; % N_rev2  
lb(6) = 3; % q1   
lb(7) = 3; % q2  
lb(8) = 0;% v_inf 
lb(9) = -pi;% alpha  
lb(10) = -pi;% beta  
lb(11) = 0;% v_infga
lb(12) = pi/2;% alphagam
lb(13) = -pi/2;% alphagap
lb(14) = pi/2;% betagam
lb(15) = -pi/2;% betagap
%upper boundary 
ub(1) = date2mjd2000([2027 1 1 0 0 0]);  % t0  
ub(2) = 300; % TOF12 
ub(3) = 800; % TOF23
ub(4) = 3; % N_rev1
ub(5) = 3;  % N_rev2
ub(6) = 7; % q1
ub(7) = 7; % q2 
ub(8) = sqrt(11);% v_inf 
ub(9) = pi;% alpha 
ub(10) = pi;% beta  
ub(11) = 4;% v_infga
ub(12) = 3/2*pi;% alphagam
ub(13) = pi/2;% alphagap
ub(14) = 3/2*pi;% betagam
ub(15) = pi/2;% betagap


Bound = [lb; ub];

options = optimoptions('gamultiobj', 'Display', 'Iter', ...
                       'PopulationSize', 100, 'StallGenLimit', 200, ... %          
                       'MaxGenerations', 20, ...
                       'ParetoFraction', 0.35, ...
                       'UseParallel', true, 'PopInitRange',Bound);
[SOL,feval,exitflag] = gamultiobj(@(x) conway_flyby(x,data_stacks), 15,[],[],[],[],lb,ub,options);

%% plots (to be fixed)
chosen = paretoplot(SOL, feval);

t0              = SOL(chosen,1); 
TOF12           = SOL(chosen,2); 
TOF23           = SOL(chosen,3);
N_rev1          = round(SOL(chosen,4));  
N_rev2          = round(SOL(chosen,5));
q1              = SOL(chosen,6);     
q2              = SOL(chosen,7);
%departure
v_inf           = SOL(chosen,8);
alpha           = SOL(chosen,9); 
beta            = SOL(chosen,10);
%gravity assist
v_infga         = SOL(chosen,11);
alphagam        = SOL(chosen,12); 
alphagap        = SOL(chosen,13); 
betagam         = SOL(chosen,14); 
betagap         = SOL(chosen,15); 

%fly by position
[kep2, muS]   = uplanet(t0 + TOF12, 3); %flyby at Earth
[R2, v2] = kep2car2(kep2, muS);  %km....

muE = 3.986004330000000e+05;
rE = 6.371010000000000e+03;

%gravity assist check
xversga = v2/norm(v2);
zversga = cross(R2, v2)/norm(cross(R2,v2));
yversga = cross(xversga, zversga);

V_infm = (cos(betagam)* (cos(alphagam)*xversga + sin(alphagam)*yversga) + sin(betagam)*zversga);
V_infp = (cos(betagap)* (cos(alphagap)*xversga + sin(alphagap)*yversga) + sin(betagap)*zversga);
delta = acos(dot(V_infm, V_infp));
e = 1/ sin(delta/2);
rp = muE/norm(v_infga)^2 *(e-1);

if rp < 1.3*rE
        printf('Unfeasible solution')
else
    
    delta_vga = v_infga*(V_infp - V_infm);
    v2 = v2 + delta_vga;

    [kep1, ~] = uplanet(t0              ,   3);
    [kep3, ~]   = uplanet(t0 + TOF12 + TOF23,   4);

    [R1, v1] = kep2car2(kep1, muS); %km....
    [R3, v3] = kep2car2(kep3, muS);  %km....

    %definition of plane of motion
    r1norm = norm(R1);
    r2norm = norm(R2);
    r3norm = norm(R3);

    r1vers = R1/r1norm;
    r2vers = R2/r2norm;
    r3vers = R3/r3norm;

    RIvcRFv12 = cross(r1vers, r2vers);
    RCRRv12 = RIvcRFv12/norm(RIvcRFv12);

    if RIvcRFv12(3) <0
        RCRRv12 = -RCRRv12;
    end

    RIvcRFv23 = cross(r2vers, r3vers);
    RCRRv23 = RIvcRFv23/norm(RIvcRFv23);

    if RIvcRFv23(3) <0
        RCRRv23 = -RCRRv23;
    end
  [ m23, T23, r23, z23, s23, vr23, vt23, vz2, acc_inplane23, acc_out23, acc23, TH23, L23, gamma123, gamma223, gamma23, v1perp23, v2perp23, v1tra23, v2tra23, vnorm23, dmdt23, T_inplane23, T_outplane23, time23, TOFr23] = ...
  Conway(TOF23, N_rev2, q2, r2norm, r3norm, r2vers, r3vers, RCRRv23, RIvcRFv23, v2, v3, muS, data_stacks);

    data_stacks.Mdry = m23(end);

    %TMI maneuver
    V_infdep = v_inf*(sin(beta)*cos(alpha)*r1vers + ...
                     sin(beta)*sin(alpha)*cross(RCRRv12,r1vers) + ...
                     cos(beta)*RCRRv12);
    v1 = v1 + V_infdep;

    [ m12, T12, r12, z12, s12, vr12, vt12, vz12, acc_inplane12, acc_out12, acc12, TH12, L12, gamma112, gamma212, gamma12, v1perp12, v2perp12, v1tra12, v2tra12, vnorm12, dmdt12, T_inplane12, T_outplane12, time12, TOFr12] = ...
        Conway(TOF12, N_rev1, q1, r1norm, r2norm, r1vers, r2vers, RCRRv12, RIvcRFv12, v1, v2, muS, data_stacks);

end
  % plots    
    
RE = zeros(length(TOFr12) + length(TOFr23), 3);
RM = RE;
REnorm = zeros(length(TOFr12) + length(TOFr23),1); 
RMnorm = REnorm;

TOFr = [TOFr12 TOFr23]; r=[r12 r23]; z=[z12 z23]; vt=[vt12 vt23];  vr=[vr12 vr23]; TH = [TH12 TH23];
rr = RE; vv = rr;

        vers = r1vers;
        RCRRv = RCRRv12;
for i =1:length(TOFr)
    kepE = uplanet(t0+TOFr(i), 3);
    kepM = uplanet(t0+TOFr(i), 4);
    RE(i,:) = kep2car2(kepE, muS);
    RM(i,:) = kep2car2(kepM, muS);
    REnorm(i) = norm(RE(i,:));
    RMnorm(i) = norm(RM(i,:));
    if i > length(TOFr12)
        vers = r2vers;
        RCRRv = RCRRv23;
    end
    [rr(i,:), vv(i,:)] = refplane2car( r(i), z(i),  vt(i), vr(i), TH(i), vers, RCRRv);
end

% 
% dirt1 = cross(RCRRv12 , r1vers);
% dirt3 = cross(RCRRv23 , r3vers);
% 
% vr1 = r1vers * vr12(1);
% vr3 = r3vers * vr23(end);
% 
% vt1 = dirt1 * vt12(1);
% vt3 = dirt3 * vt23(end);
% 
% vdep = vr1 + vt1 + v1perp*RCRRv12; varr = vr23 + vt23 + v2perp*RCRRv23;
% v_inf1 = vdep - vE0; v_inf2 = vMend - varr;
% 
% % optimization result
% datedep=mjd20002date(t0);
% datega=mjd20002date(t0 + TOF1 );
% datearr=mjd20002date(t0 + TOF1 + TOF2);
% fprintf('OPTIMIZATION RESULTS \n')
% fprintf('Departure Date \t %d %d %d \n', datedep(3), datedep(2),datedep(1));
% fprintf('Fly-by Date \t %d %d %d \n', datega(3), datega(2),datega(1));
% fprintf('Arrival Date \t %d %d %d \n', datearr(3), datearr(2),datearr(1));
% fprintf('%%%%%% Second Arc \n')
% fprintf('TOF1 \t \t %d days \n', TOF1);
% fprintf('N_rev1 \t \t %d \n', N_rev1);
% fprintf('q1 \t \t %d \n', q1);
% fprintf('%%%%%% Second Arc \n')
% fprintf('TOF2 \t \t %d days \n', TOF2);
% fprintf('N_rev2 \t \t %d \n', N_rev2);
% fprintf('q2 \t \t %d \n', q2);
% % fprintf('Departure v_inf: \t %d km/s (C3 = %d km^2/s^2) \n', norm(v_inf1), norm(v_inf1)^2)
% % fprintf('Arrival v_inf: \t %d km/s \n', norm(v_inf2))
% fprintf('Mass ratio: \t %d \n', m2(end)/m1(1))
% fprintf('Fuel Mass Fraction: \t %d \n', (m1(1) - m2(end))/m1(1))
% fprintf('Propellant mass: \t %d kg \n', m1(1) - m2(end))
% fprintf('%%%%%% FLY-BY \n')
% fprintf('rp: \t %d km \n',rp)

figure()
subplot(2,2,[1 3])
plot3(RM(:,1), RM(:,2), RM(:,3)), hold on,
plot3(RE(:,1), RE(:,2), RE(:,3)), hold on,
plot3(rr(:,1), rr(:,2), rr(:,3),'--'), hold on,
axis equal,  title('complete path')
%%
subplot(2,2,2), 
plot(TOFr, RMnorm), hold on, plot(TOFr, REnorm), hold on, plot(TOFr, r), 
hold off, title('in-plane motion')
subplot(2,2,4), 
plot(TOFr, RM*RCRRv), hold on, plot(TOFr, RE*RCRRv), hold on, plot(TOFr, z), 
hold off, title('out-of-plane motion')

figure()
sgtitle('Thrust Profile')
subplot(5,2,1), plot(TOFr, acc_inplane), title('a_{inplane}')
subplot(5,2,2), plot(TOFr, acc_out), title('a_{outplane}')
subplot(5,2,3), plot(TOFr, T_inplane), title('T_{inplane}')
subplot(5,2,4), plot(TOFr, T_outplane), title('T_{outplane}')
subplot(5,2,[5 6]), plot(TOFr, acc), title('a_{tot}')
subplot(5,2,[7 8]), plot(TOFr, T), title('T')
subplot(5,2,[9 10]), plot(TOFr, m), title('m')

figure()
sgtitle('Flight path angle')
plot(TOFr, gamma,'DisplayName','spacecraft'), hold on 
yline(gamma1,'DisplayName','Departure'); hold on, yline(gamma2,'DisplayName','Mars');
legend(),

figure()
plot(TOFr, vt12,'DisplayName','$v_{t}$'), hold on
yline(norm(v1tra),'DisplayName','$v^E_\theta$'); hold on, yline(norm(v2tra),'DisplayName','$v^M_\theta$');
plot(TOFr, v1norm,'DisplayName','$|v|$'), hold on, plot(TOFr, vr, 'DisplayName','$v_{r}$')
legend()



function chosen = paretoplot(SOL, feval)
figure()
sgtitle('Pareto Front')
minTOF = 1e5; chosen = 0;
for i = 1:length(feval)
    if feval(i,1) <= 0.22
        plot(feval(i,1), feval(i,2), 'ro','HandleVisibility','off'), hold on
    else 
        plot(feval(i,1), feval(i,2), 'k+','HandleVisibility','off'), hold on
    end
          if SOL(i,2) < minTOF
            minTOF = SOL(i,2) + SOL(i,3);
            chosen = i;
         end
end
plot(feval(chosen,1), feval(chosen,2),'go', 'DisplayName',strcat('Min TOF (', num2str(SOL(chosen,2)),')')); hold off
xlabel('max(abs(T))'), ylabel('m_P') , legend()
end