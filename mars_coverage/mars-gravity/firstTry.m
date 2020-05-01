clear, clc, close all

mu = astroConstants(14);
R = 3389.5;

SMA = 13500;
INC = deg2rad(55);

walker = [21 3 2];
bw = 20;
lon = [-170 170];
lat = [-80 80];
disc = [8 8];
alt = 0;

tic
[cov, N_mesh, N_mesh_mean] = getMinCoverage(walker, bw, SMA, INC, lon, lat, disc, alt);
toc
%%
tic

n_lon = 8;
n_lat = 8;
LON = linspace(-170, 170, n_lon);
LAT = linspace(-80, 80, n_lat);
figure
imagesc(LON, LAT, N_mesh)
title('Min sats cov - with 20 deg beamwidth', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlim([-170,170])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
% colormap(gray)


figure
imagesc(LON, LAT, N_mesh_mean)
title('Mean sats cov - 20 deg beamwidth', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlim([-170,170])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
% colormap(gray)

toc