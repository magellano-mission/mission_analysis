%{
script for analyze ECS visibility

%}

%% Matlab Initialization
clear; close all; clc

%% Figure Initialization
set(0,'DefaultFigureUnits', 'normalized');
set(0,'DefaultFigurePosition',[0 0 1 1]);
set(0,'DefaultTextFontSize',18);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultAxesXGrid','on')
set(0,'DefaultAxesYGrid','on')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');

%% 
I = imread('mars.jpg');                            % earth image
RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];                       % earth image x sizes
RI.YWorldLimits = [-90 90];                         % earth image y sizes

DateInit = [2025, 1, 1, 0, 0, 0];                          % date of start
DayInit = date2mjd2000(DateInit);                           % initial time  [days]
DateEnd = [2025, 1, 5, 0, 0, 0];                            % date of end
DayEnd = date2mjd2000(DateEnd);                             % end time  [days]
N = 10000;
days = linspace(DayInit, DayEnd, N);
dTdays = (DayEnd - DayInit);
dT = dTdays*86400;                              % simulation time [s]
t = linspace(0, dT, N);


%% ECS orbit
OrbPar(1) = 6400;
OrbPar(2) = 0;
OrbPar(3) = 0*pi/180;
OrbPar(4)  = 0;
OrbPar(5) = 0;
OrbPar1 = OrbPar; OrbPar2 = OrbPar;
theta0_sat1 = 0; theta0_sat2 = 120*pi/180;

Mars = 4;
Earth = 3;

V = 1.6318e11;                          % mars volume [km^3]
R = nthroot(3*V/(4*pi), 3);             % mars equivalent radius [km]
mi = 42828.3;                           % mars gravity constant [km^2/s^3]
n = sqrt(mi/OrbPar(1)^3);               % mars angular velocity [rad/s]
Rsun = 6.955e5;                         % Sun equivalent radius [km]
norm_rSc = OrbPar(1);

phi1 = acos(R/norm_rSc);                

los = zeros(N, 1);
RR1 = zeros(N, 3); RR2 = zeros(N, 3);
RR3 = zeros(N, 3); RR4 = zeros(N, 3); 
for i = 1:N
    OrbPar1(6) = theta0_sat1 + n*t(i);
    OrbPar2(6) = theta0_sat2 + n*t(i);
    mjd2000 = days(i);
    
    [OrbParM, ~] = uplanet(mjd2000, Mars);
    [OrbParE, ~] = uplanet(mjd2000, Earth);
    
    rS2M = kep2car_r_only(OrbParM);
    rS2E = kep2car_r_only(OrbParE); norm_rS2E = norm(rS2E);
    rSc1 = kep2car_r_only(OrbPar1);
    rSc2 = kep2car_r_only(OrbPar2);
    
    rM2E = rS2E - rS2M;      norm_rM2E = norm(rM2E);
    rS2Sc1 = rS2M + rSc1;    norm_rS2Sc1 = norm(rS2Sc1);
    rS2Sc2 = rS2M + rSc2;    norm_rS2Sc2 = norm(rS2Sc2);
    
    RR1(i, :) = rS2E;
    RR2(i, :) = rS2M;
    RR3(i, :) = rM2E/norm(rM2E)*2*OrbPar1(1);
    RR4(i, :) = rSc1;
    
    % Mars Eclipse
    phi_sat1 = acos(dot(rM2E, rSc1)/(norm_rSc*norm_rM2E));
    phi_sat2 = acos(dot(rM2E, rSc2)/(norm_rSc*norm_rM2E));
    phi2 = acos(R/norm_rM2E);
    
    % Sun Eclipse
    theta_sat1 = acos(dot(rS2E, rS2Sc1)/(norm_rS2E*norm_rS2Sc1));
    theta_sat2 = acos(dot(rS2E, rS2Sc2)/(norm_rS2E*norm_rS2Sc2));
    theta1_sat1 = acos(Rsun/norm_rS2Sc1);
    theta1_sat2 = acos(Rsun/norm_rS2Sc2);
    theta2 = acos(Rsun/norm_rS2E);
    
    if phi_sat1 < (phi1 + phi2)% || phi_sat2 < (phi1 + phi2) 
        los_Mars = 1;
    else 
        los_Mars = 0;
    end    
    
%     if theta_sat1 < (theta1_sat1 + theta2) || theta_sat2 < (theta1_sat2 + theta2)
%         los_Sun = 1;
%     else
%         los_Sun = 0;
%     end
    
    los(i) = los_Mars;%*los_Sun;
    
end

figure; plot(t/60, los, 'LineWidth', 2);
xlabel('t [days]'); ylabel('los'); ylim([-0.1, 1.1]);
title('line of sight ECS-Earth with both satellites');

% ll = (los_Mars == 1);
% nl = not(ll);
% 
%% mars-earth orbits
% figure; hold on;
% plot3(RR1(ll,1),RR1(ll,2),RR1(ll,3),'b'); plot3(RR1(nl,1),RR1(nl,2),RR1(nl,3),'r')
% plot3(RR2(ll,1),RR2(ll,2),RR2(ll,3),'b'); plot3(RR2(nl,1),RR2(nl,2),RR2(nl,3),'r')

% v = VideoWriter('Sun_POV', 'MPEG-4');
% v.FrameRate = 60;
% open(v);

% figure('units', 'normalized', 'outerposition', [0 0 1 1], 'Name', 'Sun_POV', 'NumberTitle', 'off', 'Color', 'black');
% hold on
% view(172, 2)
% set(gca,'Color','none','visible','off')
% 
% plot3(RR1(:, 1),RR1(:, 2),RR1(:, 3),'b'); 
% hold on;
% plot3(RR2(:, 1),RR2(:, 2),RR2(:, 3),'r');
% plot3(0, 0, 0, 'yo','MarkerSize', 10, 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'y')
% 
% axis equal;
% 
% h1 = plot3(RR1(1, 1), RR1(1, 2), RR1(1, 3), 'bo','MarkerSize', 5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
% h2 = plot3(RR2(1, 1), RR2(1, 2), RR2(1, 3), 'ro','MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
% 
% for i = 1:N
%     
%     h1.XData = RR1(i, 1);
%     h1.YData = RR1(i, 2);
%     h1.ZData = RR1(i, 3);
%     
%     h2.XData = RR2(i, 1);
%     h2.YData = RR2(i, 2);
%     h2.ZData = RR2(i, 3);
%     
%     drawnow
%     
% %     frame = getframe(gcf);
% %     writeVideo(v, frame);
%     
% end
% close(v);
% close('Sun_POV')

%% on mars pov video
% 
% figure('units', 'normalized', 'outerposition', [0.3 0.3 0.3 0.7], 'Name', 'Mars_POV', 'NumberTitle', 'off', 'Color', 'black');
% hold on
% view(-140, 30)
% set(gca,'Color','none','visible','off')
% set(gca,'YDir','normal')
% axis([-1.3e4, 1.3e4, -1.3e4, 1.3e4, -1.3e4, 1.3e4])
% 
% v = VideoWriter('Mars_POV');
% v.FrameRate = 60;
% v.Quality = 100;
% open(v);
% 
% [X, Y, Z] = ellipsoid(0, 0, 0, R, R, R, 100); % spheric centered earth
% planet = surf(X, Y, -Z,'Edgecolor', 'none');
% set(planet,'FaceColor','texturemap','Cdata',I)
% 
% h4 = plot3(RR4(:, 1), RR4(:, 2), RR4(:, 3), 'LineWidth', 1);
% h1 = plot3(RR4(1, 1), RR4(1, 2), RR4(1, 3), 'bo','MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
% h2 = quiver3(0, 0, 0, RR3(1, 1), RR3(1, 2), RR3(1, 3),'Color','r');
% h3 = plot3(0, 0, 0, 'ko','MarkerSize', 1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
% 
% Nv = 15*v.FrameRate;               % total frames
% 
% RR4v = RR4(1:Nv, :);
% RR3v = RR3(1:floor(N/Nv):N, :);
% 
% for i = 1:Nv
%     h1.XData = RR4v(i, 1);
%     h1.YData = RR4v(i, 2);
%     h1.ZData = RR4v(i, 3);
%     
%     hold on
%     delete(h2)
%     h2 = quiver3(0, 0, 0, RR3v(i, 1), RR3v(i, 2), RR3v(i, 3),'Color','r');
%     
%     tof = strcat('time: ', " ", num2str(floor(dTdays*i*11/N)), ' [days]');
%     legend(h3, {tof}, 'TextColor', 'w', 'Location', 'northeast')
%     
%     drawnow 
% 
%     frame = getframe(gcf);
%     writeVideo(v, frame);
%     
% end
% 
% close(v)
% close('Mars_POV')
