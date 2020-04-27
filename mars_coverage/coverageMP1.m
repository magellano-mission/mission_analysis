clear, clc, close all

mu = astroConstants(14);
R = 3389.5;

SMA = 6500;
INC = deg2rad(111);


n = 2000;
n_sat = 5;
n_orbits = 1;

n_lon = 8;
n_lat = 8;

alt = 0;

LON = linspace(-180, 180, n_lon);
LAT = linspace(-80, 80, n_lat);

T = 2 * pi * sqrt(SMA^3 / mu);
tspan = linspace(0, 3*T, n);

X0 = [SMA, 1e-8, INC, 0, 0, 0];

gamma = deg2rad(20);

opt = odeset('AbsTol', 1e-6, 'RelTol', 1e-6);

N = zeros(n,n_sat);
N_mesh = zeros(n_lat, n_lon);
N_mesh_mean = zeros(n_lat, n_lon);

tic
for lo = 1 : n_lon
    for la = 1 : n_lat
        
        for kk = 1 : n_sat
            X0(6) = (kk - 1) * pi /2;
            [T,Y] = trajectory(X0, tspan, mu, opt);
            theta = footPrintRadius(gamma, Y, alt);
            for jj = 1 : n
                N(jj, kk) = coverageNumber(LAT(la), LON(lo), T(jj), Y(jj,:), theta(jj));
            end
        end
        
        NN = sum(N,2);
        N_min = min(NN);
        N_mesh(la, lo) = N_min;
        N_mean = mean(NN);
        N_mesh_mean(la, lo) = N_mean;
    end
end
toc
%%

figure
imagesc(LON, LAT, N_mesh_mean)
title('Minimum satellites coverage on surface with 40 deg beamwidth', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
% colormap(gray)
