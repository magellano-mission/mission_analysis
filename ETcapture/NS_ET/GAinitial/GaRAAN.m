clc; clear

orbit.mi = 42828;
orbit.T = -0.2;

orbit.Om0 = 6.2832;
orbit.iF = 45*pi/180;
orbit.DOm = 240*pi/180;

orbit.a = 563.7788e+003;
orbit.e = 0.93;
orbit.om = 277.6965e-003;
orbit.m0 = 1.5965e+003;


% Output: SOL --> th0, i0, Dt
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
                       'PopulationSize', 300, 'StallGenLimit', 300, ... %          
                       'MaxGenerations', 500, ...
                       'UseParallel', true, 'PopInitRange', Bound);
                   
[SOL, feval, exitflag] = ga(@(x) BCforCapturing(x, orbit), 3, [], [], [], [], lb, ub, [], options);