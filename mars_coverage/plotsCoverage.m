close all

load('MP1 - MEAN.mat')

% figure
% imagesc(LON, LAT, N_mesh)
% title('MPI - Minimum sat - 40 deg', 'interpreter', 'latex', 'FontSize', 20)
% hold on
% yline(0)
% xlim([-180,180])
% ylim([-80, 80])
% xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
% ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
% colorbar
% % colormap(gray)
% ax = gca;
% ax.YDir = 'normal';

figure
imagesc(LON, LAT, N_mesh_mean)
title('MPI - Avarage satellites coverage', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
ax = gca;
ax.YDir = 'normal';
% colormap(gray)


load('MP2 - MEAN.mat')

% figure
% imagesc(LON, LAT, N_mesh)
% title('MPII - Avarage satellites coverage', 'interpreter', 'latex', 'FontSize', 20)
% hold on
% yline(0)
% xlim([-180,180])
% ylim([-80, 80])
% xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
% ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
% colorbar
% % colormap(gray)
% ax = gca;
% ax.YDir = 'normal';

% figure
% imagesc(LON, LAT, N_mesh_mean)
% title('MPII - Mean sat - 40 deg', 'interpreter', 'latex', 'FontSize', 20)
% hold on
% yline(0)
% xlim([-180,180])
% ylim([-80, 80])
% xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
% ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
% colorbar
% % colormap(gray)
% ax = gca;
% ax.YDir = 'normal';

load('MP2 - NO MS - MEAN.mat')

figure
imagesc(LON, LAT, N_mesh)
title('MPII - Minimum satellites coverage', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
% colormap(gray)
ax = gca;
ax.YDir = 'normal';

figure
imagesc(LON, LAT, N_mesh_mean)
title('MPII - Avarage satellites coverage',  'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
% colormap(gray)
ax = gca;
ax.YDir = 'normal';

load('MP2 - MS - 9sat - MEAN.mat')

% figure
% imagesc(LON, LAT, N_mesh)
% title('MPII - 9 sat - Minimum sat - 40 deg', 'interpreter', 'latex', 'FontSize', 20)
% hold on
% yline(0)
% xlim([-180,180])
% ylim([-80, 80])
% xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
% ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
% colorbar
% % colormap(gray)
% ax = gca;
% ax.YDir = 'normal';

% figure
% imagesc(LON, LAT, N_mesh_mean)
% title('MPII - 9 sat - Mean sat - 40 deg', 'interpreter', 'latex', 'FontSize', 20)
% hold on
% yline(0)
% xlim([-180,180])
% ylim([-80, 80])
% xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
% ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
% colorbar
% % colormap(gray)
% % ax = gca;
% % ax.YDir = 'normal';

load('MP3 - MEAN.mat')

figure
imagesc(LON, LAT, N_mesh)
title('MPIII - Minimum sat - 40 deg', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
% colormap(gray)
ax = gca;
ax.YDir = 'normal';

% figure
% imagesc(LON, LAT, N_mesh_mean)
% title('MPIII - Avarage satellites coverage', 'interpreter', 'latex', 'FontSize', 20)
% hold on
% yline(0)
% xlim([-180,180])
% ylim([-80, 80])
% xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
% ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
% colorbar
% % colormap(gray)
% ax = gca;
% ax.YDir = 'normal';

load('MP3 - NO MS - MEAN.mat')

figure
imagesc(LON, LAT, N_mesh)
title('MPIII - Minimum satellites coverage', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
% colormap(gray)
ax = gca;
ax.YDir = 'normal';

figure
imagesc(LON, LAT, N_mesh_mean)
title('MPIII - Avarage satellites coverage', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
% colormap(gray)
ax = gca;
ax.YDir = 'normal';

load('MP3 - NO MS - 1000km - MEAN.mat')

% figure
% imagesc(LON, LAT, N_mesh)
% title('MPIII - NO MS - 1000km - Minimum sat - 40 deg', 'interpreter', 'latex', 'FontSize', 20)
% hold on
% yline(0)
% xlim([-180,180])
% ylim([-80, 80])
% xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
% ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
% colorbar
% % colormap(gray)
% ax = gca;
% ax.YDir = 'normal';

figure
imagesc(LON, LAT, N_mesh_mean)
title('MPIII at 1000km - Avarage satellites coverage', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
% colormap(gray)
ax = gca;
ax.YDir = 'normal';

% Senza mars stationary - fatto 
% Meno di 12 sat (9 ?) - fatto con e senza mars stationary
% MPI
% Diverso angolo - da fare 
% Rifare la configurazione MP2 a 500km. (fatti tutti eccetto MP2)
% Fare solo MPII a 1000km.
% MP3 - fatto

% TODO:
%  MP1 - done
%  MP3 FULL 
%  MP3 FULL 1000km
%  CAMBRIARE ANGOLO - done

%% 

figure
sgtitle('Avarage satellites coverage', 'interpreter', 'latex')

load('MP1 - MEAN.mat')
subplot(1,3,1)
imagesc(LON, LAT, N_mesh_mean)
title('MPI', 'interpreter', 'latex', 'FontSize', 12)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 12)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 12)
colorbar
ax = gca;
ax.YDir = 'normal';
% colormap(gray)

load('MP2 - NO MS - MEAN.mat')
subplot(1,3,2)
imagesc(LON, LAT, N_mesh_mean)
title('MPII',  'interpreter', 'latex', 'FontSize', 12)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 12)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 12)
colorbar
% colormap(gray)
ax = gca;
ax.YDir = 'normal';

load('MP3 - NO MS - MEAN.mat')
subplot(1,3,3)
imagesc(LON, LAT, N_mesh_mean)
title('MPIII', 'interpreter', 'latex', 'FontSize', 12)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 12)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 12)
colorbar
% colormap(gray)
ax = gca;
ax.YDir = 'normal';