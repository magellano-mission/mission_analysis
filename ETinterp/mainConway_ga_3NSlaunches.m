%Conway-based approach
run('MatlabSetup.m')
%% SPACECRAFT CHARACTERISTICS
run('Stacks_config.m')
dataNS = dataNS1;

choice = menu('Do you want to run a new simulation or load a SOL?','Run a simulation','Load SOL');
if choice == 1 || choice == 0
%% BOUNDARIES
% t0, TOF, N_rev, q, v_inf, alpha, beta, v_infcap, alphacap, betacap 
%%%%%
% lower boundary
lb = zeros(1,21); ub = lb;
lb(1)   = date2mjd2000([2024 1 1 0 0 0]); %t0
lb(2)   = 700; %TOF
lb(3)   = 0; %N_rev
lb(4)   = 3;% q
lb(5)   = 0;% v_inf dep
lb(6)   = 0;% alpha dep
lb(7)   = 0;% beta dep
lb(8)   = 700; %TOF
lb(9)   = 0; %N_rev
lb(10)  = 3;% q
lb(11)  = 0;% v_inf dep
lb(12)  = 0;% alpha dep
lb(13)  = 0;% beta dep
lb(14)  = 700; %TOF
lb(15)  = 0; %N_rev
lb(16)  = 3;% q
lb(17)  = 0;% v_inf dep
lb(18)  = 0;% alpha dep
lb(19)  = 0;% beta dep

% upper boundary 
ub(1)   = date2mjd2000([2026 1 1 0 0 0]);
ub(2)   = 1000;
ub(3)   = 3;
ub(4)   = 8;
ub(5)   = 3; %v_inf dep
ub(6)   = pi;% alpha dep
ub(7)   = pi;% beta dep
ub(8)   = 1000;
ub(9)   = 3;
ub(10)  = 8;
ub(11)  = 3; %v_inf dep
ub(12)  = pi;% alpha dep
ub(13)  = pi;% beta dep
ub(14)  = 1100;
ub(15)  = 3;
ub(16)  = 8;
ub(17)  = 3; %v_inf dep
ub(18)  = pi;% alpha dep
ub(19)  = pi;% beta dep

%flexibility margins
lb(20)  = -5;
lb(21)  = -5;

ub(20)  = 5;
ub(21)  = 5;
 
Bound = [lb; ub];

%% choose optimization method
optimization_method = menu('Remember to set boudaries... \n optimization method:','ga','gamultiboj','particleswarm');

switch optimization_method
    case 1
    options = optimoptions('ga', 'Display', 'Iter', ...
                           'PopulationSize', 200, 'StallGenLimit', 200, ...         
                           'MaxGenerations', 200, ...
                           'UseParallel', true, 'PopInitRange',Bound);
    [SOL,feval,exitflag] = ga(@(x) ga_conway_3NSlaunches(x,dataNS), 21,[],[],[],[],lb,ub,[], [3 9 15],options);
    chosen = 1;
    case 2
    options = optimoptions('gamultiobj', 'Display', 'Iter', ...
                           'PopulationSize', 200, 'StallGenLimit', 200, ...         
                           'MaxGenerations', 200, ...
                           'ParetoFraction', 0.35, ...
                           'UseParallel', true, 'PopInitRange',Bound);
    [SOL,feval,exitflag] = gamultiobj(@(x) gamultiobj_conway_3NSlaunches(x,dataNS), 21,[],[],[],[],lb,ub,[],options);
    [chosen, chosenT] = paretoplot(SOL, feval,dataNS); %min TOF chosen
    chosen = chosenT;
    case 3
    options = optimoptions('particleswarm', 'SwarmSize', 1000, 'MaxIterations', 20, 'Display', 'Iter','UseParallel', true);
    [SOL,feval,exitflag] = particleswarm(@(x) ga_conway_3NSlaunches(x, dataNS), 21, lb, ub, options);
    chosen = 1;
    case 0
        fprintf('\n loaded a solution... \n')
    load('SOLpartswarm2.mat')
    chosen = 1;
end
else
    load('SOLpartswarm2.mat')
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

%% first launch
diary results
fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%%%')
fprintf('\n')

run ('launch1_3NS.m') 

diary off
fprintf('\n press any button to continue... \n')
pause()
%% second launch
diary on

run('launch2_3NS.m')

diary off
fprintf('\n press any button to continue... \n')
pause()
%% third launch
diary on

run('launch3_3NS.m')

fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%%%')
fprintf('\n')
diary off


%% LOCAL UTILITY FUNCTIONS

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
