function plotCoverage(Y, t, theta)

[xx,yy,zz,~,~,~] = cone(rad2deg(theta),[0 0 0; Y(1) Y(2) Y(3)], norm(Y));

x = xx(:,end);
y = yy(:,end);
z = zz(:,end);

r = [x, y, z];
tt = ones(length(x), 1) * t;

[lat, lon] = groundTrack(r, tt);

lon = real(lon);

hold on
plot(lon, lat, '.b')


end