
clc;clear;close all; 
set(0, 'DefaultTextInterpreter', 'Latex');
set(0, 'DefaultLegendInterpreter', 'Latex');
%% Relative Translational Dynamics Simulation 
par = struct;
par.a  = 7170*1e3;
par.e  = 0.05;
par.mu = 3.98e14;
par.i  = deg2rad(15); 
par.OM = deg2rad(0); 
par.om = deg2rad(340);
par.JT = [500;550;600];
par.JL = [500;550;600];

T = 2*pi*sqrt(par.a^3/par.mu);
rho0    = [25;25;50]; 
rhoDot0 = [0.0;-0.0555;0.0];
th0     = deg2rad(20) ;
thDot0  = deg2rad(0.0656);
w0     = deg2rad([0.006;0.006;0.12]);
q0      = [0;0;0;1];
P10     = [1.5;1.5;0];
P20     = [-1.5;-1.5;0];


X0 = [rho0;rhoDot0;th0];
par.opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,X] = ode113(@relativeTranslation, [0 T], X0, par.opt, par);
[~,XCW] = ode113(@ClohessyWiltshire, t, X0(1:end-1), par.opt, par);
figure()
subplot(3,1,1)
plot(t,[X(:,1),XCW(:,1)], 'LineWidth',1.4);grid on;xlim([0 t(end)])
xlabel('t [s]'); ylabel ('x [m]')
legend('$\rho_{CM}$')
ax = gca; ax.FontSize = 12; 
subplot(3,1,2)
plot(t,X(:,2), 'LineWidth',1.4);grid on;xlim([0 t(end)])
xlabel('t [s]'); ylabel ('y [m]')
ax = gca; ax.FontSize = 12; 
subplot(3,1,3)
plot(t,X(:,3), 'LineWidth',1.4); grid on;xlim([0 t(end)])
xlabel('t [s]'); ylabel ('z [m]')
ax = gca; ax.FontSize = 12; 
sgtitle('Relative Translational Dynamics')

%% Coupled Relative Rotational and Translational Dynamics 
X02 =[rho0;rhoDot0;w0;q0;th0; P10;P20];
[t,X2] = ode113(@coupledRD,[0 T], X02, par.opt, par);
rhoCM = X2(:,1:3);

rho1 = zeros(length(X2),3);
rho2 = zeros(length(X2),3);

parout2 = zeros(length(X2), length(X02(15:end)));
for kk = 1 : length(X2)
    [~, parout2(kk,:)] = coupledRD(t(kk), X2(kk,:)', par);
    rho1 (kk,:) = parout2(kk,1:3);
    rho2 (kk,:) = parout2(kk,4:6);

end

figure()
subplot(3,1,1)
plot(t,rhoCM(:,1), 'LineWidth',1.2);grid on;xlim([0 t(end)]); hold on ; 
plot(t,rho1(:,1),'r','LineWidth',1.2);
plot(t,rho2(:,1),'k','LineWidth',1);
xlabel('t [s]'); ylabel ('x [m]')
ax = gca; ax.FontSize = 12; 
legend('$\rho_{CM}$','$\rho_{10}$','$\rho_{20}$', 'FontSize',16)
subplot(3,1,2)
plot(t,rhoCM(:,2), 'LineWidth',1.2);grid on;xlim([0 t(end)]);hold on;
plot(t,rho1(:,2),'r','LineWidth',1.2);
plot(t,rho2(:,2),'k','LineWidth',1);
xlabel('t [s]'); ylabel ('y [m]')
ax = gca; ax.FontSize = 12; 
subplot(3,1,3)
plot(t,rhoCM(:,3), 'LineWidth',1.2); grid on;xlim([0 t(end)]);hold on;
plot(t,rho1(:,3),'r','LineWidth',1.2);
plot(t,rho2(:,3),'k','LineWidth',1);
xlabel('t [s]'); ylabel ('z [m]')
ax = gca; ax.FontSize = 12; 
sgtitle('Relative Rotational and Translational Dynamics')


figure()
subplot(3,1,1)
plot(t,rho1(:,1)-rhoCM(:,1), 'LineWidth',1.4);grid on;xlim([0 t(end)]); hold on ; 
plot(t,rho2(:,1)-rhoCM(:,1), 'r', 'LineWidth',1.4);
xlabel('t [s]'); ylabel ('x [m]')
legend('$\Delta\rho_{1}$','$\Delta\rho_{2}$')
ax = gca; ax.FontSize = 12; 
subplot(3,1,2)
plot(t,rho1(:,2)-rhoCM(:,2), 'LineWidth',1.4);grid on;xlim([0 t(end)]); hold on ; 
plot(t,rho2(:,2)-rhoCM(:,2), 'r', 'LineWidth',1.4);
xlabel('t [s]'); ylabel ('y [m]')
ax = gca; ax.FontSize = 12; 
subplot(3,1,3)
plot(t,rho1(:,3)-rhoCM(:,3), 'LineWidth',1.4);grid on;xlim([0 t(end)]); hold on ; 
plot(t,rho2(:,3)-rhoCM(:,3), 'r', 'LineWidth',1.4);
xlabel('t [s]'); ylabel ('z [m]')
ax = gca; ax.FontSize = 12; 
sgtitle('Feature Points Dynamics')
