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
 
data_stacks.Isp = 5000;                                      % specific impulse [s]
data_stacks.Mdry = 2800;                                      % Total Mass of the s/c [kg]
data_stacks.n_int = 1000;
%%%%%
% t0, TOF, N_rev, q, v_inf, alpha, beta, v_infcap, alphacap, betacap 
%%%%%
%lower boundary
lb = zeros(1,5); ub = lb;
lb(1) = date2mjd2000([2024 1 1 0 0 0]); %t0
lb(2) = 600; %TOF
lb(3) = 0; %N_rev
lb(4) = 0;% v_inf dep
lb(5) = 0;% alpha dep
% lb(6) = 0;% v_inf arr
% lb(7) = 0;% alpha arr
% lb(8) = 0;% beta arr
%upper boundary 
ub(1) = date2mjd2000([2027 1 1 0 0 0]);
ub(2) = 1400;
ub(3) = 3;
ub(4) = sqrt(9); %v_inf dep
ub(5) = pi;% alpha dep
% ub(6) = 0.1;% v_inf arr
% ub(7) = pi;% alhpa arr
% ub(8) = pi;% beta arr
 
Bound = [lb; ub];
% load('optim1.m')
%optimization of thrust only
options = optimoptions('ga', 'Display', 'Iter', ...
                       'PopulationSize', 200, 'StallGenLimit', 10, ... %          
                       'MaxGenerations', 200, ...
                       'UseParallel', true, 'PopInitRange',Bound);% 'InitialPopulationMatrix', optim1);
[SOL,feval,exitflag] = ga(@(x) ga_conway2D(x,data_stacks), 10,[],[],[],[],lb,ub,[],[],options);

% options = optimoptions('gamultiobj', 'Display', 'Iter', ...
%                        'PopulationSize', 200, 'StallGenLimit', 200, ... %          
%                        'MaxGenerations', 200, 'ParetoFraction', 0.35, ...
%                        'UseParallel', true, 'PopInitRange',Bound);
% [SOL,feval,exitflag] = gamultiobj(@(x) ga_conway2D(x,data_stacks), 10,[],[],[],[],lb,ub,options);
% feval = feval/1e4;
%% plots
% chosen = paretoplot(SOL, feval)
% chosen = length(feval)
chosen = find(feval(:,1)==min(feval(:,1)))
% chosen = 1;
data_stacks.n_int = 1000;
t0         =          SOL(chosen,1); 
TOF        =          SOL(chosen,2);
N_rev      =    round(SOL(chosen,3));
v_inf      =          SOL(chosen,4);
alpha      =          SOL(chosen,5);

[kepEarth, muS] = uplanet(t0      ,3);
[kepMars, ~]    = uplanet(t0 + TOF,4);
 
[R1, v1] = kep2car2(kepEarth, muS); %km....
[R2, v2] = kep2car2(kepMars, muS);  %km....
 
%definition of plane of motion
r1norm = norm(R1);
r2norm = norm(R2);
 
r1vers = R1/r1norm;
r2vers = R2/r2norm;
 
hh = cross(r1vers, r2vers);
hvers = hh/norm(hh);
 
if hh(3) <0
    hvers = -hvers;
end

%adding TMI maneuver
v1 = v1 + v_inf*(cos(alpha)*r1vers + sin(alpha)*cross(hvers, r1vers));
    
v1 = v1 + v_inf1;
             
% %adding v_inf at Mars capture
% v2 = v2 + v_infcap*(sin(betacap)*cos(alphacap)*r2vers + ...
%                  sin(betacap)*sin(alphacap)*cross(RCRRv,r2vers) + ...
%                  cos(betacap)*RCRRv);             
             
 
[ m, T, r, vr, vt, acc_inplane, acc, TH, L, gamma1, gamma2, gamma, v1perp, v2perp, v1tra, v2tra, vnorm, dmdt, T_inplane, time, TOFr] = ...
    Conway2D(TOF, N_rev, r1norm, r2norm, r1vers, r2vers, hvers, hh, v1, v2, muS, data_stacks);
 
if ~isnan(r) 
    
    
[kepEarth, muS] = uplanet(t0, 3);
[rE0, vE0] = kep2car2(kepEarth, muS);
[kepMars, muS] = uplanet(t0 + TOF, 4);
[rMend, vMend] = kep2car2(kepMars, muS);
 
r1vers = rE0/norm(rE0);
r2vers = rMend/norm(rMend) ; 
    
R1 = zeros(length(TOFr), 3); RM = R1;
REnorm = zeros(length(TOFr),1); RMnorm = REnorm;
rr = R1; vv = rr;
for i =1:length(TOFr)
    kepE = uplanet(t0+TOFr(i), 3);
    kepM = uplanet(t0+TOFr(i), 4);
    R1(i,:) = kep2car2(kepE, muS);
    RM(i,:) = kep2car2(kepM, muS);
    REnorm(i) = norm(R1(i,:));
    RMnorm(i) = norm(RM(i,:));
    [rr(i,:), vv(i,:)] = refplane2car( r(i), 0,  vt(i), vr(i), TH(i), r1vers, hvers);
end
 
dirt1 = cross(hvers , r1vers);
dirt2 = cross(hvers , r2vers);
 
vr1 = r1vers * vr(1);
vr2 = r2vers * vr(end);
 
vt1 = dirt1 * vt(1);
vt2 = dirt2 * vt(end);
 
v1 = vr1 + vt1 + v_inf1; v2 = vr2 + vt2;
v_inf2 = vMend - v2;
 
% optimization result
datedep=mjd20002date(t0);
datearr=mjd20002date(t0 + TOF);
fprintf('OPTIMIZATION RESULTS \n')
fprintf('Departure Date \t %d %d %d \n', datedep(3), datedep(2),datedep(1));
fprintf('Arrival Date \t %d %d %d \n', datearr(3), datearr(2),datearr(1));
fprintf('TOF \t \t %d days (%d yrs)\n', TOF, TOF/365);
fprintf('N_rev \t \t %d \n', N_rev);
fprintf('Departure v_inf: \t %d km/s (C3 = %d km^2/s^2) \n', norm(v_inf1), norm(v_inf1)^2)
fprintf('Arrival v_inf: \t %d km/s (C3 = %d km^2/s^2) \n', norm(v_inf2), norm(v_inf2)^2)
fprintf('Mass ratio: \t %d \n', m(end)/m(1))
fprintf('Fuel Mass Fraction: \t %d \n', (m(1) - m(end))/m(1))
fprintf('Propellant mass: \t %d kg \n', m(1) - m(end))
fprintf('Max T: \t %d N \n',feval(chosen,1)) 
 
figure()
subplot(2,1,1)
plot3(RM(:,1), RM(:,2), RM(:,3)), hold on,
plot3(R1(:,1), R1(:,2), R1(:,3)), hold on,
plot3(rr(:,1), rr(:,2), rr(:,3),'--'), hold on,
axis equal,  title('complete path')

subplot(2,1,2), 
sgtitle('Conway 2D')
plot(TOFr, RMnorm), hold on, plot(TOFr, REnorm), hold on, plot(TOFr, r), 
hold off, title('in-plane motion')
 
figure()
sgtitle('Thrust Profile - planar trajectory')
subplot(3,2,1), plot(TOFr, acc_inplane), title('a_{inplane}')
subplot(3,2,3), plot(TOFr, T_inplane), title('T_{inplane}')
subplot(3,2,2), plot(TOFr, acc), title('a_{tot}')
subplot(3,2,4), plot(TOFr, T), title('T')
subplot(3,2,[5 6]), plot(TOFr, m), title('m')

figure()
sgtitle('Flight path angle')
plot(TOFr, gamma,'DisplayName','spacecraft'), hold on 
yline(gamma1,'DisplayName','Departure'); hold on, yline(gamma2,'DisplayName','Mars');
legend(),

figure()
plot(TOFr, vt,'DisplayName','$v_{t}$'), hold on
yline(norm(v1tra),'DisplayName','$v^E_\theta$'); hold on, yline(norm(v2tra),'DisplayName','$v^M_\theta$');
plot(TOFr, vnorm,'DisplayName','$|v|$'), hold on, plot(TOFr, vr, 'DisplayName','$v_{r}$')
legend()

else
    fprintf('No real solution for Conway algorithm \n')
end


function chosen = paretoplot(SOL, feval)
figure()
sgtitle('Pareto Front')
minTOF = 1e5; chosen = 0;
threshold = 0.25; %[N]
for i = 1:length(feval)
    if feval(i,1) <= threshold
        plot(feval(i,1), feval(i,2), 'ro','HandleVisibility','off'), hold on
    else 
        plot(feval(i,1), feval(i,2), 'k+','HandleVisibility','off'), hold on
    end
         if SOL(i,2) < minTOF
            minTOF = SOL(i,2);
            chosen = i;
         end
end
plot(feval(chosen,1), feval(chosen,2),'go', 'DisplayName',strcat('Min TOF (', num2str(SOL(chosen,2)),')')); hold off
xline(threshold,'r'), xlabel('max(abs(T))'), ylabel('m_P') , legend()
end