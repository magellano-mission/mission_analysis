clc; clear
orbit.Om0 = 0;
orbit.iF = 45*pi/180;
orbit.DOm = 120*pi/180;
orbit.a = 5.2939e+06;
orbit.e = 0.99;
orbit.om = 0.03154;
orbit.T = -0.2;
orbit.m0 = 2399.4;
orbit.mi = 42828;

% lower boundary
lb = zeros(1, 3); ub = lb;

lb(1) = 0;
lb(2) = 0.1*pi/180;
lb(3) = 1e4;
% upper boundary 
ub(1) = 2*pi;
ub(2) = 179*pi/280;
ub(3) = 86400*80;

Bound = [lb; ub];

% options = optimoptions('gamultiobj', 'Display', 'Iter', ...
options = optimoptions('ga', 'Display', 'Iter', ...
                       'PopulationSize', 100, 'StallGenLimit', 200, ... %          
                       'MaxGenerations', 200, ...
                       'UseParallel', true, 'PopInitRange', Bound);
                   
[SOL, feval, exitflag] = ga(@(x) BCforCapturing(x, orbit), 3, [], [], [], [], lb, ub, [], options);