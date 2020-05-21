%Conway-based approach

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

%definition of interplanetary arc
TOF = 3.5*365;%[days]
N_rev = 2; %max 2-3 for Conway hp on convergence
q = 3; %qmin = 3

data_stacks.Isp = 4300;                                      % specific impulse [s]
data_stacks.Mdry = 8000;                                      % Total Mass of the s/c [kg]

%lower boundary
lb(1) = date2mjd2000([2023 1 1 0 0 0]);
lb(2) = 300;
lb(3) = 0;
lb(4) = 3;
%upper boundary 
ub(1) = date2mjd2000([2023+3 1 1 0 0 0]);
ub(2) = 1500;
ub(3) = 3;
ub(4) = 10;

options = optimoptions('gamultiobj', 'Display', 'Iter', 'MaxGenerations', 1000, 'UseParallel', true);
[SOL,feval,exitflag] = gamultiobj(@(x) ga_conway(x,data_stacks), 4,[],[],[],[],lb,ub,options);
t0 = SOL(1); 
TOF = SOL(2);
N_rev = SOL(3);
q = SOL(4);

%capture (computation of v_arrival necessary)
%%%%%%

%plots
if ~isnan(r) 
RE = zeros(length(TOFr), 3); RM = RE;
REnorm = zeros(length(TOFr),1); RMnorm = REnorm;
rr = RE; vv = rr;
for i =1:length(TOFr)
    kepE = uplanet(data_stacks.t0sym+TOFr(i), 3);
    kepM = uplanet(data_stacks.t0sym+TOFr(i), 4);
    RE(i,:) = kep2car2(kepE, muS);
    RM(i,:) = kep2car2(kepM, muS);
    REnorm(i) = norm(RE(i,:));
    RMnorm(i) = norm(RM(i,:));
    [rr(i,:), vv(i,:)] = refplane2car( r(i), z(i), rE0/norm(rE0), vr(i), vt(i), RCRRv, TH(i), muS); %da controllare
end

r1vers = rE0/norm(rE0);
r2vers = rMend/norm(rMend) ;
dirt1 = cross(RCRRv , r1vers);
dirt2 = cross(RCRRv , r2vers);

vr1 = r1vers * vr(1);
vr2 = r2vers * vr(end);

vt1 = dirt1 * vt(1);
vt2 = dirt2 * vt(end);

v1 = vr1 + vt1 + v1perp; v2 = vr2 + vt2 + v2perp;
v_inf1 = v1 - vE0; v_inf2 = vMend - v2;

fprintf('Departure v_inf: \t %d \n', norm(v_inf1))
fprintf('Arrival v_inf: \t %d \n', norm(v_inf2))

figure()

plot3(r.*cos(TH), r.*sin(TH), z), hold on, axis equal
plot3(REnorm(1), 0, 0, 'o'), hold on,
plot3(cos(L)*RMnorm(end), RMnorm(end)*sin(L), 0, 'o'), hold on,


figure()
sgtitle('Thrust Profile')
subplot(3,1,1), plot(TOFr, acc), title('a_{tot}')
subplot(3,1,2), plot(TOFr, T), title('T')
subplot(3,1,3), plot(TOFr, m), title('mass')

figure()
sgtitle('Flight path angle')
plot(TOFr, gamma), hold on 
yline(gamma1); hold on, yline(gamma2);

figure()
plot(TOFr, vt,'DisplayName','$v_{t}$'), hold on
yline(norm(v1tra),'DisplayName','$v^E_\theta$'); hold on, yline(norm(v2tra),'DisplayName','$v^M_\theta$');
plot(TOFr, vnorm,'DisplayName','$|v|$'), hold on, plot(TOFr, vr, 'DisplayName','$v_{r}$')
legend()
else
    fprintf('No real solution for Conway algorithm \n')
end

figure()
plot(TOFr, z)
