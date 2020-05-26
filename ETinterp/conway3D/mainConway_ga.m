%Conway-based approach
 
%definition of parameters
close all, clear, clc
% Figure Initialization    
run('PlotsInitialization.m')

data_stacks.Isp = 4300;                                       % specific impulse [s]
data_stacks.Mdry = 2400;                                      % Total Mass of the s/c [kg]
data_stacks.n_int = 1000;
%%%%%
% t0, TOF, N_rev, q, v_inf, alpha, beta, v_infcap, alphacap, betacap 
%%%%%
%lower boundary
lb = zeros(1,7); ub = lb;
lb(1) = date2mjd2000([2024 1 1 0 0 0]); %t0
lb(2) = 600; %TOF
lb(3) = 1; %N_rev
lb(4) = 3;% q
lb(5) = 0;% v_inf dep
lb(6) = 0;% alpha dep
lb(7) = 0;% beta dep

%upper boundary 
ub(1) = date2mjd2000([2027 1 1 0 0 0]);
ub(2) = 1400;
ub(3) = 3;
ub(4) = 8;
ub(5) = sqrt(9); %v_inf dep
ub(6) = pi;% alpha dep
ub(7) = pi;% beta dep

 
Bound = [lb; ub];
% load('optim1.m')
%optimization of thrust only
% options = optimoptions('ga', 'Display', 'Iter', ...
%                        'PopulationSize', 200, 'StallGenLimit', 200, ... %          
%                        'MaxGenerations', 200, ...
%                        'UseParallel', true, 'PopInitRange',Bound);% 'InitialPopulationMatrix', optim1);
% [SOL,feval,exitflag] = ga(@(x) ga_conway(x,data_stacks), 10,[],[],[],[],lb,ub,[],[],options);

% %%multiobj optimization of thrust and mass
options = optimoptions('gamultiobj', 'Display', 'Iter', ...
                       'PopulationSize', 200, 'StallGenLimit', 200, ... %          
                       'MaxGenerations', 200, ...
                       'ParetoFraction', 0.35, ...
                       'UseParallel', true, 'PopInitRange',Bound);
[SOL,feval,exitflag] = gamultiobj(@(x) gamultiobj_conway(x,data_stacks), 7,[],[],[],[],lb,ub,[],options);
feval(:,1) = feval(:,1)/10000;
%% solution choice 
chosen = paretoplot(SOL, feval); %min TOF chosen
% chosen = length(feval) %last value chosen
chosen = find(feval(:,1)==min(feval(:,1))) %min T chosen
% chosen = 13
data_stacks.n_int = 1000;
t0         =          SOL(chosen,1); 
TOF        =          SOL(chosen,2);
N_rev      =    round(SOL(chosen,3));
q          =          SOL(chosen,4);
v_inf      =          SOL(chosen,5);
alpha      =          SOL(chosen,6);
beta       =          SOL(chosen,7);

%% propagation to check that the computed states are compliant
[kepEarth, muS] = uplanet(t0      ,3);
[kepMars, ~]    = uplanet(t0 + TOF,4);
 
[R1, v1] = kep2car2(kepEarth, muS); %km....
[R2, v2] = kep2car2(kepMars, muS);  %km....
 
%definition of plane of motion
r1norm = norm(R1);
r2norm = norm(R2);
 
r1vers = R1/r1norm;
r2vers = R2/r2norm;
 
RIvcRFv = cross(r1vers, r2vers);
RCRRv = RIvcRFv/norm(RIvcRFv);
 
if RIvcRFv(3) <0
    RCRRv = -RCRRv;
end
 
%adding TMI maneuver
v1 = v1 + v_inf*(sin(beta)*cos(alpha)*r1vers + ...
                 sin(beta)*sin(alpha)*cross(RCRRv,r1vers) + ...
                 cos(beta)*RCRRv);          
             
 
[ m, T, r, z, s, vr, vt, vz, acc_inplane, acc_out, acc, TH, L, gamma1, gamma2, gamma, v1perp, v2perp, v1tra, v2tra, vnorm, dmdt, T_inplane, T_outplane, th_dot, time, TOFr] = ...
    Conway(TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, RCRRv, RIvcRFv, v1, v2, muS, data_stacks);



if ~isnan(r) 
    X0 = [r(1), TH(1), z(1), vr(1), th_dot(1), vz(1), m(1)];
    X = zeros(length(TH), 7);
    X(1,:) = X0;
    %RK4 forward integration
    for i = 1:length(TH)-1
        x = X(i,:);
        
        hhh = time(i+1) - time(i);
        K1 = propagate_conway(time(i), x, acc_inplane(i), acc_out(i), gamma(i), muS, data_stacks);
        K2 = propagate_conway(time(i)+ hhh/2, x + hhh/2*K1, 0.5*(acc_inplane(i+1) + acc_inplane(i)), 0.5*(acc_out(i+1) + acc_out(i)), 0.5*(gamma(i+1)-gamma(i)), muS, data_stacks);
        K3 = propagate_conway(time(i)+ hhh/2, x + hhh/2*K2, 0.5*(acc_inplane(i+1) + acc_inplane(i)), 0.5*(acc_out(i+1) + acc_out(i)), 0.5*(gamma(i+1)-gamma(i)), muS, data_stacks);
        K4 = propagate_conway(time(i)+ hhh, x + hhh*K3, acc_inplane(i+1), acc_out(i+1), gamma(i+1), muS, data_stacks);    
        
        X(i+1,:) = X(i,:) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end
    
TX_inplane = acc_inplane .*X(:,7)' * 1000;
TX_outplane = acc_out .*X(:,7)' * 1000;
TX = (TX_inplane.^2 + TX_outplane.^2).^0.5;
figure()
sgtitle('states validation')
subplot(8,2,1), plot(TOFr, r,':'), subplot(8,2,2),plot(TOFr, X(:,1),':'), hold on
subplot(8,2,3), plot(TOFr, TH, ':'), subplot(8,2,4), plot(TOFr, X(:,2),':'), hold on
subplot(8,2,5), plot(TOFr, z, ':'), subplot(8,2,6), plot(TOFr, X(:,3),':'), hold on
subplot(8,2,7), plot(TOFr, vr, ':'), subplot(8,2,8), plot(TOFr, X(:,4), ':'), hold on
subplot(8,2,9), plot(TOFr, th_dot, ':'), subplot(8,2,10), plot(TOFr, X(:,5), ':'), hold on
subplot(8,2,11), plot(TOFr, vz, ':'), subplot(8,2,12), plot(TOFr, X(:,6), ':'), hold on
subplot(8,2,13), plot(TOFr, m, ':'), subplot(8,2,14), plot(TOFr, X(:,7), ':'), hold on
subplot(8,2,15), plot(TOFr, T, ':'), subplot(8,2,16), plot(TOFr, TX, ':'), hold on

else
    fprintf('No real solution coming from Conway algorithm \n')
end
%% plots

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
    [rr(i,:), vv(i,:)] = refplane2car( r(i), z(i),  vt(i), vr(i), TH(i), r1vers, RCRRv);
end
 
 
dirt1 = cross(RCRRv , r1vers);
dirt2 = cross(RCRRv , r2vers);
 
vr1 = r1vers * vr(1);
vr2 = r2vers * vr(end);
 
vt1 = dirt1 * vt(1);
vt2 = dirt2 * vt(end);
 
v1 = vr1 + vt1 + v1perp*RCRRv; v2 = vr2 + vt2 + v2perp*RCRRv;
v_inf1 = v1 - vE0; v_inf2 = vMend - v2;
 
% optimization result
datedep=mjd20002date(t0);
datearr=mjd20002date(t0 + TOF);
fprintf('OPTIMIZATION RESULTS \n')
fprintf('Departure Date \t %d %d %d \n', datedep(3), datedep(2),datedep(1));
fprintf('Arrival Date \t %d %d %d \n', datearr(3), datearr(2),datearr(1));
fprintf('TOF \t \t %d days (%d yrs)\n', TOF, TOF/365);
fprintf('N_rev \t \t %d \n', N_rev);
fprintf('q \t \t %d \n \n', q);
fprintf('Departure v_inf: \t %d km/s (C3 = %d km^2/s^2) \n', norm(v_inf1), norm(v_inf1)^2)
fprintf('Arrival v_inf: \t %d km/s (C3 = %d km^2/s^2)   \n \n', norm(v_inf2), norm(v_inf2)^2)
fprintf('Isp: \t %d s\n', data_stacks.Isp)
fprintf('Mdry: \t %d kg\n', data_stacks.Mdry)
fprintf('Propellant mass: \t %d kg \n', m(1) - m(end))
fprintf('Mass ratio: \t %d \n', m(end)/m(1))
fprintf('Fuel Mass Fraction: \t %d \n', (m(1) - m(end))/m(1))
fprintf('max T : \t %d N\n',max(abs(T)))
figure()
sgtitle('')
subplot(2,2,[1 3])
plot3(RM(:,1), RM(:,2), RM(:,3)), hold on,
plot3(R1(:,1), R1(:,2), R1(:,3)), hold on,
plot3(rr(:,1), rr(:,2), rr(:,3),'--'), hold on,
axis equal,  title('complete path')
 
subplot(2,2,2), 
plot(TOFr, RMnorm), hold on, plot(TOFr, REnorm), hold on, plot(TOFr, r), 
hold off, title('in-plane motion')
subplot(2,2,4), 
plot(TOFr, RM*RCRRv), hold on, plot(TOFr, R1*RCRRv), hold on, plot(TOFr, z), 
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
plot(TOFr, vt,'DisplayName','$v_{t}$'), hold on
yline(norm(v1tra),'DisplayName','$v^E_\theta$'); hold on, yline(norm(v2tra),'DisplayName','$v^M_\theta$');
plot(TOFr, vnorm,'DisplayName','$|v|$'), hold on, plot(TOFr, vr, 'DisplayName','$v_{r}$')
legend()
%  annotation(gcf,'textbox',...
%     [0.15 0.75 0.29 0.15],...
%     'String',{strcat('Dep :',num2str(datedep(3)),'.', num2str(datedep(2)),'.', num2str(datedep(1)), '...', ...
%                      'Arr :', num2str(datearr(3)),'.', num2str(datearr(2)),'.', num2str(datearr(1))), ...
%               strcat('Nrev :',num2str(N_rev),'...','TOF :',num2str(TOF),'days'), ...
%               strcat('v_{inf} dep :', num2str(v_inf1),'km/s', 'v_{inf} arr :',num2str(v_inf2), 'km/s'), ...
%               strcat('max|T| :', num2str(max(abs(T))),'N','...','Isp :',num2str(data_stacks.Isp),'s','...','Mdry :', num2str(data_stacks.Mdry),'kg')},...
%     'FitBoxToText','on');
else
    fprintf('No real solution for Conway algorithm \n')
end
 

%%
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
