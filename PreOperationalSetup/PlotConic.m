function [h, R] = PlotConic(OrbPar, Color, LineWidth, LineType, th)

N = length(th);
R = zeros(N, 3);
for ii = 1:N
    OrbPar(6) = th(ii)*pi/180;
    R(ii, :) = kep2car_r_only(OrbPar);
end
hold on;
h = plot3(R(:, 1), R(:, 2), R(:, 3), LineType, 'Color', Color, 'LineWidth', LineWidth );