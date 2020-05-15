function [h, R] = FakeRAANplot(OrbPar, Color, LineWidth, LineType, th, RAAN)

N = length(th);
R = zeros(N, 3);
for ii = 1:N
    OrbPar(4) = RAAN(ii);
    OrbPar(6) = th(ii)*pi/180;
    R(ii, :) = kep2car_r_only(OrbPar);
end
    
h = plot3(R(:, 1), R(:, 2), R(:, 3), LineType, 'Color', Color, 'LineWidth', LineWidth );