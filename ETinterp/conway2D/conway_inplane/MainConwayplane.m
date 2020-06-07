%%
clear; close all; clc

data_stacks.Isp = 6000;                                      % specific impulse [s]
data_stacks.Mdry = 2400;                                      % Total Mass of the s/c [kg]
data_stacks.n_int = 1000;

% lower boundary
lb = zeros(1, 5); ub = lb;
lb(1) = date2mjd2000([2024 1 1 0 0 0]);
lb(2) = 700;
lb(3) = 0;
lb(4) = 0;
lb(5) = 0;
% upper boundary 
ub(1) = date2mjd2000([2025 1 1 0 0 0]);
ub(2) = 1500;
ub(3) = 3;
ub(4) = 2;
ub(5) = pi/4;

Bound = [lb; ub];

options = optimoptions('ga', 'Display', 'Iter', ...
                       'PopulationSize', 100, 'StallGenLimit', 200, ... %          
                       'MaxGenerations', 20, ...
                       'UseParallel', true, 'PopInitRange', Bound);
                   
[SOL, feval, exitflag] = ga(@(x) ConwayOptim(x, data_stacks), 5, [], [], [], [], lb, ub, [], 3, options);

%% plot solution
data_stacks.n_int = 1000;
t0         =          SOL(1); 
TOF        =          SOL(2);
N_rev      =    round(SOL(3));
v_inf      =          SOL(4);
alpha      =          SOL(5);

beta = pi/2;
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
v_inf1 = v_inf*(sin(beta)*cos(alpha)*r1vers + ...
                 sin(beta)*sin(alpha)*cross(hvers,r1vers) + ...
                 cos(beta)*hvers);
v1 = v1 + v_inf1;
                     
 
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
fprintf('Max T: \t %d N \n',feval(1)) 
 
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

[y, m, T, r, TH, Vinf] = ConwayOptim(SOL, data_stacks);

figure()
subplot(2,1,1), plot(m)
subplot(2,1,2), plot(T)

figure()
plot(cos(TH).*r, sin(TH).*r), axis equal
else
    fprintf('No real solution for Conway algorithm \n')
end

