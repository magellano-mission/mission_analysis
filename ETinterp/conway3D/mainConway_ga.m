%Conway-based approach
run('MatlabSetup.m')
%% SPACECRAFT CHARACTERISTICS
data_stacks.Isp = 4300;                                       % specific impulse [s]
data_stacks.Mdry = 1600;                                      % Total Mass of the s/c [kg]
data_stacks.n_int = 1000;
%% BOUNDARIES
% t0, TOF, N_rev, q, v_inf, alpha, beta, v_infcap, alphacap, betacap 
%%%%%
%lower boundary
% lb = zeros(1,7); ub = lb;
% lb(1) = date2mjd2000([2024 1 1 0 0 0]); %t0
% lb(2) = 600; %TOF
% lb(3) = 0; %N_rev
% lb(4) = 3;% q
% lb(5) = 0;% v_inf dep
% lb(6) = 0;% alpha dep
% lb(7) = 0;% beta dep
%upper boundary 
% ub(1) = date2mjd2000([2027 1 1 0 0 0]);
% ub(2) = 1000;
% ub(3) = 2;
% ub(4) = 8;
% ub(5) = sqrt(9); %v_inf dep
% ub(6) = pi;% alpha dep
% ub(7) = pi;% beta dep

%lower boundary
lb = zeros(1,6); ub = lb;
lb(1) = 600; %TOF
lb(2) = 0; %N_rev
lb(3) = 3;% q
lb(4) = 0;% v_inf dep
lb(5) = 0;% alpha dep
lb(6) = 0;% beta dep
%upper boundary 
ub(1) = 1000;
ub(2) = 2;
ub(3) = 8;
ub(4) = sqrt(9); %v_inf dep
ub(5) = pi;% alpha dep
ub(6) = pi;% beta dep

%lower boundary
% lb(1) = 0; lb(2) = 3;
%upper boundary
% ub(1) = 3; ub(2) = 8;
 
Bound = [lb; ub];
% load('optim1.m')
%optimization of thrust only
% options = optimoptions('ga', 'Display', 'Iter', ...
%                        'PopulationSize', 200, 'StallGenLimit', 200, ... %          
%                        'MaxGenerations', 200, ...
%                        'UseParallel', true, 'PopInitRange',Bound);% 'InitialPopulationMatrix', optim1);
% [SOL,feval,exitflag] = ga(@(x) ga_conway(x,data_stacks), 10,[],[],[],[],lb,ub,[],[],options);
%% GENETIC ALGORITHM
% %%multiobj optimization of thrust and mass
options = optimoptions('gamultiobj', 'Display', 'Iter', ...
                       'PopulationSize', 200, 'StallGenLimit', 200, ... %          
                       'MaxGenerations', 200, ...
                       'ParetoFraction', 0.35, ...
                       'UseParallel', true, 'PopInitRange',Bound);
[SOL,feval,exitflag] = gamultiobj(@(x) gamultiobj_conway(x,data_stacks), 6,[],[],[],[],lb,ub,[],options);
feval(:,1) = feval(:,1)/10000;
%% SOLUTIONS
data_stacks.Tmax = 0.2; % N
chosen = paretoplot(SOL, feval,data_stacks); %min TOF chosen
% chosen = length(feval) %last value chosen
chosen = find(feval(:,1)==min(feval(:,1))) %min T chosen
% chosen = find(feval(:,2)==min(feval(:,2))) %min mP chosen
% chosen = 1
%%
data_stacks.n_int = 1000;
%results
% t0         =          SOL(chosen,1); 
% TOF        =          SOL(chosen,2);
% N_rev      =    round(SOL(chosen,3));
% q          =          SOL(chosen,4);
% v_inf      =          SOL(chosen,5);
% alpha      =          SOL(chosen,6);
% beta       =          SOL(chosen,7);

% %results NS2
% TOF        =          SOL(chosen,1);
% N_rev      =    round(SOL(chosen,2));
% q          =          SOL(chosen,3);
% v_inf      =          SOL(chosen,4);
% alpha      =          SOL(chosen,5);
% beta       =          SOL(chosen,6);
% tf = 9.619974167516721e+03 + 8.361779091278747e+02 + 180; 
% t0 = tf - TOF;

% %results NS3
% TOF        =          SOL(chosen,1);
% N_rev      =    round(SOL(chosen,2));
% q          =          SOL(chosen,3);
% v_inf      =          SOL(chosen,4);
% alpha      =          SOL(chosen,5);
% beta       =          SOL(chosen,6);
% tf = 9.619974167516721e+03 + 8.361779091278747e+02 + 360; 
% t0 = tf - TOF;

%results NS (presentation)
v_inf      =          1.019932213771344;
alpha      =          1.339209562214939;
beta       =          2.689364771171173;
t0 = 9.618881860211146e+03; TOF = 9.084518115422542e+02;
q = 5.140866; N_rev = 1;

% %% starting from NS departure
% v_inf      =          1.019932213771344;
% alpha      =          1.339209562214939;
% beta       =          2.689364771171173;
% 
% N_rev      =    round(SOL(chosen,1));
% q          =          SOL(chosen,2);
% % N_rev = 1;
% % q = 6.652955;
% t0 = 9.618881860211146e+03;
%TOF = 9.084518115422542e+02;
% % 

%% propagation to check that the computed states are compliant

[kepEarth, muS] = uplanet(t0,3);
[kepMars, ~]    = uplanet(t0 + TOF,4);
 
[R1, v1] = kep2car2(kepEarth, muS); %km....
[R2, v2] = kep2car2(kepMars, muS);  %km....
 
%definition of plane of motion
r1norm = norm(R1);
r2norm = norm(R2);
 
r1vers = R1/r1norm;
r2vers = R2/r2norm;
 
hh = cross(r1vers, r2vers);
href = hh/norm(hh);
 
if hh(3) <0
    href = -href;
end
 
%adding TMI maneuver
v1 = v1 + v_inf*(sin(beta)*cos(alpha)*r1vers + ...
                 sin(beta)*sin(alpha)*cross(href,r1vers) + ...
                 cos(beta)*href);          
             
[ m, T, r, z, s, vr, vt, vz, acc_inplane, acc_out, acc, TH, L, gamma1, gamma2, gamma, v1perp, v2perp, v1tra, v2tra, vnorm, dmdt, T_inplane, T_outplane, th_dot, time, TOFr] = ...
    Conway(TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, href, hh, v1, v2, muS, data_stacks);

if ~isnan(r) 
    X0 = [r(1), TH(1), z(1), vr(1), th_dot(1), vz(1), m(1)];
    X = zeros(length(TH), 7);
    X(1,:) = X0;
    %RK4 forward integration
    for k = 1:length(TH)-1
        x = X(k,:);
        
        hhh = time(k+1) - time(k);
        K1 = propagate_conway(time(k), x, acc_inplane(k), acc_out(k), muS, data_stacks);
        K2 = propagate_conway(time(k)+ hhh/2, x + hhh/2*K1, 0.5*(acc_inplane(k+1) + acc_inplane(k)), 0.5*(acc_out(k+1) + acc_out(k)), muS, data_stacks);
        K3 = propagate_conway(time(k)+ hhh/2, x + hhh/2*K2, 0.5*(acc_inplane(k+1) + acc_inplane(k)), 0.5*(acc_out(k+1) + acc_out(k)), muS, data_stacks);
        K4 = propagate_conway(time(k)+ hhh, x + hhh*K3, acc_inplane(k+1), acc_out(k+1), muS, data_stacks);    
        
        X(k+1,:) = X(k,:) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end
    
TX_inplane = acc_inplane .*X(:,7)' * 1000;
TX_outplane = acc_out .*X(:,7)' * 1000;
TX = (TX_inplane.^2 + TX_outplane.^2).^0.5;
figure()
sgtitle('States validation - Dynamics states')
subplot(6,2,1), plot(TOFr, r,'k:'), ylabel('r [km]'), title('Conway Solution')
subplot(6,2,2), plot(TOFr, X(:,1),'k:'), ylabel('r [km]'), title('Propagated Solution')
subplot(6,2,3), plot(TOFr, TH, 'k:'), ylabel('$\theta$ [km]'), 
subplot(6,2,4), plot(TOFr, X(:,2),'k:'), ylabel('$\theta$ [km]')
subplot(6,2,5), plot(TOFr, z, 'k:'),  ylabel('z[km]')
subplot(6,2,6), plot(TOFr, X(:,3),'k:'), ylabel('z[km]')
subplot(6,2,7), plot(TOFr, vr, 'k:'), ylabel('$v_r$ [km]')
subplot(6,2,8), plot(TOFr, X(:,4), 'k:'), ylabel('$v_r$ [km]')
subplot(6,2,9), plot(TOFr, th_dot, 'k:'),  ylabel('$\dot{\theta}$ [km]')
subplot(6,2,10), plot(TOFr, X(:,5), 'k:'), ylabel('$\dot{\theta}$[km]')
subplot(6,2,11), plot(TOFr, vz, 'k:'),  xlabel('TOFr [days]'), ylabel('$v_z$ [km]') 
subplot(6,2,12), plot(TOFr, X(:,6), 'k:'),xlabel('TOFr [days]'), ylabel('$v_z$ [km]')

figure(), sgtitle('States validation - Mass and Thrust')
subplot(2,2,1), plot(TOFr, m, 'k:'), xlabel('TOF [days]'), ylabel('Spacecraft Mass [kg]')
title('Conway Solution'), 
subplot(2,2,2), plot(TOFr, X(:,7), 'k:'), xlabel('TOF [days]'), ylabel('Spacecraft Mass [kg]')
title('Propagated Solution')
subplot(2,2,3), plot(TOFr, T, 'k:'),  xlabel('TOF [days]'), ylabel('Thrust [N]')
subplot(2,2,4), plot(TOFr, TX, 'k:'),  xlabel('TOF [days]'), ylabel('Thrust [N]')
else
    fprintf('No real solution coming from Conway algorithm \n')
end
%% plots
run('conway3Dplots.m')
%% Direct Transcript Method
run('directTranscription.m')
%%
function chosen = paretoplot(SOL, feval, data)
figure()
sgtitle('Pareto Front')
minTOF = 1e5; chosen = 0;
for i = 1:length(feval)
    if feval(i,1) <= data.Tmax
        plot(feval(i,1), feval(i,2), 'ro','HandleVisibility','off'), hold on
    else 
        plot(feval(i,1), feval(i,2), 'k+','HandleVisibility','off'), hold on
    end
         if SOL(i,2) < minTOF
            minTOF = SOL(i,2);
            chosen = i;
         end
end
chosenT = find(feval(:,1) == min(feval(:,1)));
plot(feval(chosen,1), feval(chosen,2),'go', 'DisplayName',strcat('Min TOF (', num2str(SOL(chosen,2)),'days )')); hold on
plot(feval(chosenT,1), feval(chosenT,2),'bo', 'DisplayName',strcat('Min ( max(T(t))) (', num2str(feval(chosenT,1)),'N )')); hold off
xline(data.Tmax,'r','DisplayName','Max Thrust < 0.2 N'), xlabel('max(|T(t)|) [N]'), ylabel('m_P [kg]') , legend()
end
