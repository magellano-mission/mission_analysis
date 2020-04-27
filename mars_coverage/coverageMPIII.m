clear, clc, close all

mu = astroConstants(14);
R = 3389.5;

SMA = 11500;
INC = deg2rad(55);

P_mars = 24*3600 + 39 * 60;
SMA_stationary =  (mu * (P_mars / (2*pi))^2)^(1/3);


TT = 21;
P = 3;
F = 2;

n = 2000;
n_sat = TT/P;
n_orbits = P;

n_lon = 12;
n_lat = 12;

alt = 0;

LON = linspace(-180, 180, n_lon);
LAT = linspace(-90, 90, n_lat);


T = 2 * pi * sqrt(SMA^3 / mu);
tspan = linspace(0, 3*T, n);

X0 = [SMA, 1e-8, INC, 0, 0, 0];

gamma = deg2rad(20);

opt = odeset('AbsTol', 1e-10, 'RelTol', 1e-10);

N = zeros(n,n_sat * n_orbits + 3);
N_mesh = zeros(n_lat, n_lon);

tic
for lo = 1 : n_lon
    for la = 1 : n_lat
        
        for ii = 1 : n_orbits
            for kk = 1 : n_sat
                X0(6) = (kk - 1) * pi /2 + (ii - 1) * F * 2 * pi / TT;
                X0(4) = (ii - 1) * 2 * pi / P;
                [T,Y] = trajectory(X0, tspan, mu, opt);
                theta = footPrintRadius(gamma, Y, alt);
                for jj = 1 : n
                    N(jj, (ii-1)*n_sat + kk) = coverageNumber(LAT(la), LON(lo), T(jj), Y(jj,:), theta(jj));
                end
            end
        end
        
%         X01 = [SMA_stationary, 1e-8, 1e-15, 0, 0, 0];
%         for kk = 1 : 3
%             X01(6) = (kk - 1) * 2 * pi /3;
%             [T,Y] = trajectory(X01, tspan, mu, opt);
%             theta = footPrintRadius(gamma, Y, alt);
%             for jj = 1 : n
%                 N(jj, n_orbits * n_sat + kk) = coverageNumber(LAT(la), LON(lo), T(jj), Y(jj,:), theta(jj));
%             end
%         end
        
        NN = sum(N,2);
        N_min = min(NN);
        N_mesh(la, lo) = N_min;
    end
end
toc
%%

figure
imagesc(LON, LAT, N_mesh)
title('Minimum satellites coverage on surface with 40 deg beamwidth', 'interpreter', 'latex', 'FontSize', 20)
hold on
yline(0)
xlim([-180,180])
ylim([-80, 80])
xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
colorbar
% colormap(gray)
