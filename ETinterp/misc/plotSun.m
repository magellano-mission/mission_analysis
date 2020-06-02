function plotSun()

I = imread('Sun.jpg'); RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];  RI.YWorldLimits = [-90 90]; 
rSun = 20*almanac('Sun','Radius','kilometers','sphere');
[XSun, YSun, ZSun] = ellipsoid(0, 0, 0, rSun, rSun, rSun, 100); % spheric centered Mars
planet = surf(XSun, YSun, -ZSun,'Edgecolor', 'none','HandleVisibility','off'); hold on
set(planet,'FaceColor','texturemap','Cdata',I), axis equal
end

