% Coverage: Figures of Merit

%HP:
%circular orbits
%all constellation parameters fixed 
%
%timestep constant in all integrations
%groundtrack based approach
%conical shaped beam

%FoM: TO BE ADDED
% Percent coverage
% maximum gap
% mean gap
% time average gap
% mean response gap

% K_A = 2*pi;
% lambda0 = acos(rM/(norm(Y)));
% IAA = K_A*(1-cos(lambda0));


%% Setup
clear, close, clc

walker = [18,3,2];         % full constellation - original
% walker = [1,1,1];           % one satellite plot
% walker = [15,3,2];        % reduced constellation

bw = 40;                    % NS beamwidth [deg]
% bw = 20;                  % RS beamwidth [deg]

% Parameters
SMA = 11500;                % semi-major axis [km]
INC = deg2rad(55);          % inclination [deg]
lon = [-180, 180];          
lat = [-90,90];
disc = [50 50];             % discretization grid
alt = 0;                    % altitude at which evaluate the coverage ( ground level = 0)
timesteps = 1000;               
N_orbits = 3;               % number of orbits

%orbital periods computation
tic
[YYY, T, THETA, H] = const_orbits(walker, bw, SMA, INC, timesteps, N_orbits, alt);
toc

%% Coverage history
tic
[time_map, LON, LAT] = time_mapping(walker, YYY, T, THETA, H, lon, lat, disc);
toc
%% post-processing & figures of merit computation
trashold = 1; %number of satellites seen from grid targets
tic
[N_min, N_mean, cov, percent_cov, max_cov_gap, mean_cov_gap] = getMinCoverage(disc, time_map, T(end)-T(1), trashold);
toc

%% PLOTS
close all
figure

%mars texture reading 
texture = imread('Mars.jpg');
x=[-180 180]; y=[-90 90];
texture = flipud(texture);

colormap default
image(x, y, texture), hold on,
view(0,-90);
axis([-180 180,-90,90])
xlabel('Longitude [deg]'), ylabel('Latitude [deg]')

sss = pcolor(LON, LAT, N_min'); hold on
sss.FaceColor = 'Interp';
sss.EdgeColor = 'interp';
sss.FaceAlpha = 0.4;

% Rover positions:
plot(LON(19),LAT(32),'wo','MarkerFace','r')       % Viking 1
text(LON(19)+3,LAT(32)-3,'Viking 1','Color','w')

plot(LON(42),LAT(38),'wo','MarkerFace','r')        % Viking 2
text(LON(42)+3,LAT(38)-3,'Viking 2','Color','w')

plot(LON(25),LAT(25),'wo','MarkerFace','r')        % Opportunity
text(LON(25)+3,LAT(25)-3,'Opportunity','Color','w')

plot(LON(49),LAT(22),'wo','MarkerFace','r')        % Spirit
text(LON(49)-12,LAT(22)-3,'Spirit','Color','w')

plot(LON(44),LAT(27),'wo','MarkerFace','r')        % Insight
text(LON(44)+3,LAT(27)-3,'InSight','Color','w')

plot(LON(44),LAT(24),'wo','MarkerFace','r')        % Curiosity
text(LON(44)+3,LAT(24)-3,'Curiosity','Color','w')

plot(LON(48),LAT(5),'wo','MarkerFace','r')        % Mars Polar Lander
text(LON(48)-36,LAT(5)+3,'Mars Polar Lander','Color','w')

plot(LON(22),LAT(30),'wo','MarkerFace','r')        % ExoMars
text(LON(22)+3,LAT(30)-3,'ExoMars','Color','w')

plot(LON(3),LAT(3),'wo','MarkerFace','r')        % Skipper (Mars-Penguin)
text(LON(3)+3,LAT(3)+3,'Skipper (Mars-Penguin)','Color','w')


for i=1:length(LON)
    for j=1:length(LAT)
        plot(LON(i),LAT(j),'r+'), hold on
    end
end
title('Minimum satellites coverage on surface', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
%colorbar

%% percent_cov, max_cov_gap, mean_cov_gap
figure('Name', 'percent coverage')
image(x, y, texture), hold on,
view(0,-90);
axis([-180 180,-90,90])
xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
sss1 = pcolor(LON, LAT, percent_cov'); hold on
sss1.FaceColor = 'Interp';
sss1.EdgeColor = 'interp';
sss1.FaceAlpha = 0.7;
for i=1:length(LON)
    for j=1:length(LAT)
        plot(LON(i),LAT(j),'r+'), hold on
    end
end
title('Percent coverage on surface', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar

figure('Name','max coverage map')

image(x, y, texture), hold on,
view(0,-90);
axis([-180 180,-90,90])
xlabel('Longitude [deg]'), ylabel('Latitude [deg]')

sss = pcolor(LON, LAT, max_cov_gap'); hold on
sss.FaceColor = 'Interp';
sss.EdgeColor = 'interp';
sss.FaceAlpha = 0.7;

for i=1:length(LON)
    for j=1:length(LAT)
        plot(LON(i),LAT(j),'r+'), hold on
    end
end
title('Maximum coverage gap', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar

figure('Name','Mean coverage gap')

image(x, y, texture), hold on,
view(0,-90);
axis([-180 180,-90,90])
xlabel('Longitude [deg]'), ylabel('Latitude [deg]')

sss = pcolor(LON, LAT, mean_cov_gap'); hold on
sss.FaceColor = 'Interp';
sss.EdgeColor = 'interp';
sss.FaceAlpha = 0.7;

for i=1:length(LON)
    for j=1:length(LAT)
        plot(LON(i),LAT(j),'r+'), hold on
    end
end
title('Minimum satellites coverage on surface with 40 deg beamwidth', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
%%
figure
axis equal

image(x, y, texture), hold on,
view(0,-90);
axis([-180 180,-90,90])
xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
timevarying = pcolor(LON, LAT, time_map(:,:,1)');
colorbar
timevarying.FaceColor = 'Interp';
timevarying.EdgeColor = 'interp';
timevarying.FaceAlpha = 0.2;
% timevarying = plot(lonSS(1), latSS(1),'x');

for jjjj = 1:timesteps
%     timevarying.XData = lonSS(jjjj); hold on
%     timevarying.YData = latSS(jjjj); hold on
%     xlim([lon(1) lon(2)])
%     ylim([lat(1) lat(2)])
timevarying.CData =  time_map(:,:,jjjj)';
    drawnow
    pause(0.07)
end

%% Constellation plot & groundtrack
figure

mu = astroConstants(14);
TT = walker(1);
P = walker(2);
F = walker(3);
n_sat = TT/P;
n_orbits = P;
rM = almanac('mars','radius','kilometers','sphere');

subplot(1,2,1), 
I = imread('Mars.jpg');                            % Mars image
RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];                       % Mars image x sizes
RI.YWorldLimits = [-90 90];                         % Mars image y sizes
[X, Y, Z] = ellipsoid(0, 0, 0, rM, rM, rM, 100); % spheric centered Mars
planet = surf(X, Y, -Z,'Edgecolor', 'none');
hold on
set(planet,'FaceColor','texturemap','Cdata',I)
axis equal
set(gca,'Color','black')

subplot(1,2,2),

x=[-180 180]; y=[-90 90];
texture = flipud(texture);
image(x, y, texture), hold on,
view(0,-90);
axis([-180 180,-90,90])
xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
pointY = animatedline('Color','y','Marker','o','MarkerSize',100);
%%
for jjjj=1:timesteps
 for ii = 1 : n_orbits
     colors = ['g','y','b'];
            for kk = 1 : n_sat
                Y = squeeze(YYY((ii-1)*n_sat + kk,:,:));
                theta = THETA((ii-1)*n_sat + kk,:);
                
                subplot(1,2,1),plot3(Y(:,1),Y(:,2),Y(:,3), colors(ii)), hold on
                addpoints(pointY, Y(jjjj,1),Y(jjjj,2),Y(jjjj,3)), hold on
                
                subplot(1,2,2),
                [latSS, lonSS, ~] = groundTrack(Y, T);
                sca = scircle1(latSS(jjjj), lonSS(jjjj), rad2deg(theta(1)));
                geo = geoshow(sca(:,1),sca(:,2)); hold on
                plot(lonSS,latSS,strcat(colors(ii),'.'),'MarkerSize',0.5)
            end
 end
 pause(0.1)
  subplot(1,2,1)
% clearpoints(pointY)
end