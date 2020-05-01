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

% walker =[21,3,2];         % full constellation - original
walker = [1,1,1];           % one satellite plot
% walker = [15,3,2];        % reduced constellation

bw = 40;                    % NS beamwidth [deg]
% bw = 20;                  % RS beamwidth [deg]

% Parameters
SMA = 11500;                % semi-major axis [km]
INC = deg2rad(55);          % inclination [deg]
lon = [-180, 180];          
lat = [-90,90];
disc = [10 10];             % discretization grid
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
tic
[N_min, N_mean, cov] = getMinCoverage(disc, time_map);
toc
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

%mars texture reading 
texture = imread('Mars.jpg');
x=[-180 180]; y=[-90 90];
texture = flipud(texture);
image(x, y, texture), hold on,
view(0,-90);
axis([-180 180,-90,90])
xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
pointY = animatedline('Color','y','Marker','o','MarkerSize',100);

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
%% PLOTS
figure

x=[-180 180]; y=[-90 90];
texture = flipud(texture);
image(x, y, texture), hold on,
view(0,-90);
axis([-180 180,-90,90])
xlabel('Longitude [deg]'), ylabel('Latitude [deg]')

sss = pcolor(LON, LAT, N_min'); hold on
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

figure
image(x, y, texture), hold on,
view(0,-90);
axis([-180 180,-90,90])
xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
sss1 = pcolor(LON, LAT, N_mean'); hold on
sss1.FaceColor = 'Interp';
sss1.EdgeColor = 'interp';
sss1.FaceAlpha = 0.7;
for i=1:length(LON)
    for j=1:length(LAT)
        plot(LON(i),LAT(j),'r+'), hold on
    end
end
title('Mean satellites coverage on surface with 40 deg beamwidth', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
% xlim([-180,180])
% ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
% colormap(gray)
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