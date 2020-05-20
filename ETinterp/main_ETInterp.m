% Interplanetary travel with Electric Thrusters

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

data_stacks.t0sym = date2mjd2000([2024 1 1 0 0 0]);
data_stacks.tmax = date2mjd2000([2030 1 1 0 0 0]);

data_stacks.isPerturbed = 0;
data_stacks.isInterp = 1;
data_stacks.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
data_stacks.Isp = 10000;                                      % specific impulse [s]
data_stacks.M0 = 7000;                                      % Total Mass of the s/c [kg]
data_stacks.Across_sun = 10;                                % Cross area related to the sun [m^2]

data_stacks.event = 1;
data_stacks.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12, 'Event', @event_interp_LT);

Thrust = 0.8;                                               % value of thrust [N]
data_stacks.T = [Thrust; 0; 0];                           % thrust [N] (@TNH) (constant profile for the moment)

i_e = 0;
% while i_e ~= 2
data_stacks.t0sym = data_stacks.t0sym + 1;
% initial point
X0 = zeros(7,1);
[kepEarth, muS] = uplanet(data_stacks.t0sym,3);
[X0(1:3), X0(4:6)] = kep2car2(kepEarth, muS);
[T, Y, parout, t_e, y_e, i_e] = cart_cont_thrust_model(X0, data_stacks);
% end
TOF = (T(end) - T(1))/86400/365;
fprintf('TOF %d days \n', TOF)
rE = zeros(360,3);
rM = zeros(360,3);
rE0 = kep2car2(kepEarth, muS);
x = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options);

[kepMars, ~] = uplanet(T(end)/86400,4);
rMend = kep2car2(kepMars, muS);
for d = 1:360
    kepEarth(6) = deg2rad(d);
    kepMars(6) = deg2rad(d);
    rE(d, :) = kep2car2(kepEarth, muS);
    rM(d, :) = kep2car2(kepMars, muS);
end
figure()
plot3 (rE(:,1), rE(:,2), rE(:,3),'--'), hold on
plot3 (rM(:,1), rM(:,2), rM(:,3),'--'), hold on
plot3 (rE0(1), rE0(2), rE0(3),'o'), hold on
plot3 (rMend(1), rMend(2), rMend(3),'o'), hold on
plot3(Y(:,1), Y(:,2), Y(:,3),'k'), axis equal
%save TOF, thrust, delta v 


