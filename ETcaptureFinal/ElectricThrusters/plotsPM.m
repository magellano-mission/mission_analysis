%% Plot for report
clc; clear; close all
% NS1 NS2 NS3 ECS(42) ECS(52) RSEQ RSPOL 
X = categorical({'NS1', 'NS2', 'NS3', 'ECS 42', 'ECS 52', 'RS1 42', 'RS1 52', 'RS2 42', 'RS2 52'});
X = reordercats(X,{'NS1', 'NS2', 'NS3', 'ECS 42', 'ECS 52', 'RS1 42', 'RS1 52', 'RS2 42', 'RS2 52'});

SMs_mass = [55.3; 50.4; 46.5; 62.4; 59.4; 71; 66;...
    71.8; 59.5];

SMs_tof = [400; 369.7; 306.6; 579.4; 450; 644.3; 522.4;...
     636.5; 472.8];

figure()
subplot(2,1,1)
SMs1 = bar(X, SMs_mass);
SMs1.FaceColor = 'flat';
SMs1.CData(4,:) = [0.9490, 0.4745, 0.3137];
SMs1.CData(6,:) = [0.9490, 0.4745, 0.3137];
SMs1.CData(8,:) = [0.9490, 0.4745, 0.3137];
SMs1.CData(1,:) = [0.1020, 0.6667, 0.74120];
SMs1.CData(2,:) = [0.1020, 0.6667, 0.74120];
SMs1.CData(3,:) = [0.1020, 0.6667, 0.74120];
SMs1.CData(5,:) = [0.1020, 0.6667, 0.74120];
SMs1.CData(7,:) = [0.1020, 0.6667, 0.74120];
SMs1.CData(9,:) = [0.1020, 0.6667, 0.74120];
ylabel('$\Delta M_{prop}$ [kg]','Interpreter','Latex')


subplot(2,1,2)
SMs2 = bar(X, SMs_tof);
SMs2.FaceColor = 'flat';
SMs2.CData(4,:) = [0.9490, 0.4745, 0.3137];
SMs2.CData(6,:) = [0.9490, 0.4745, 0.3137];
SMs2.CData(8,:) = [0.9490, 0.4745, 0.3137];
SMs2.CData(1,:) = [0.1020, 0.6667, 0.74120];
SMs2.CData(2,:) = [0.1020, 0.6667, 0.74120];
SMs2.CData(3,:) = [0.1020, 0.6667, 0.74120];
SMs2.CData(5,:) = [0.1020, 0.6667, 0.74120];
SMs2.CData(7,:) = [0.1020, 0.6667, 0.74120];
SMs2.CData(9,:) = [0.1020, 0.6667, 0.74120];
ylabel('ToF [days]','Interpreter','Latex')


