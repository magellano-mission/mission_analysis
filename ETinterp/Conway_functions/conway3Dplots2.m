function conway3Dplots2(t0, TOFr, T_cart, rr, vv, muS)
%output: cartesian components of states
%conway 3d plots
   
% [kepEarth, muS] = uplanet(t0, 3);
% [rE0, vE0] = kep2car2(kepEarth, muS);
% [kepMars, muS] = uplanet(t0 + TOF, 4);
% [rMend, vMend] = kep2car2(kepMars, muS);

AU = astroConstants(2);
TU = (AU^3/muS).^0.5;
    
RE = zeros(length(TOFr), 3); RM = RE;
REnorm = zeros(length(TOFr),1); RMnorm = REnorm;

for i =1:length(TOFr)
    kepE = uplanet(t0+TOFr(i), 3);
    kepM = uplanet(t0+TOFr(i), 4);
    RE(i,:) = kep2car2(kepE, muS);
    RM(i,:) = kep2car2(kepM, muS);
    REnorm(i) = norm(RE(i,:));
    RMnorm(i) = norm(RM(i,:));
end


figure()
sgtitle('')
rmp = plot3(RM(:,1), RM(:,2), RM(:,3),'-','HandleVisibility','off'); hold on,
rep = plot3(RE(:,1), RE(:,2), RE(:,3),'-','HandleVisibility','off'); hold on,
rrp = plot3(rr(:,1), rr(:,2), rr(:,3),'--','HandleVisibility','off'); hold on,
rmpend = plot3(RM(end,1), RM(end,2), RM(end,3),'o','DisplayName','Mars @ Arrival'); hold on,
rep1 = plot3(RE(1,1), RE(1,2), RE(1,3),'o','DisplayName','Earth @ Departure'); hold on,
rrp1 = plot3(rr(1,1), rr(1,2), rr(1,3),'+','HandleVisibility','off'); hold on,
rrpend = plot3(rr(end,1), rr(end,2), rr(end,3),'+','HandleVisibility','off');hold on

%velocity plot
for i = 1:1:length(TOFr)
    quivV= quiver3(rr(i,1), rr(i,2), rr(i,3), vv(i,1), vv(i,2), vv(i,3), 2e6, 'Color','k', 'HandleVisibility','off'); hold on
    quivT = quiver3(rr(i,1), rr(i,2), rr(i,3), T_cart(i,1), T_cart(i,2), T_cart(i,3), 2e8,'Color','r', 'HandleVisibility','off'); hold on

end
quivV.HandleVisibility = 'on';
quivV.DisplayName = 'Velocity Direction';
quivT.HandleVisibility = 'on';
quivT.DisplayName = 'Thrust Direction';

plotSun()
legend()
rmpend.Color = rmp.Color; rep1.Color = rep.Color; rrp1.Color = rrp.Color; rrpend.Color = rrp.Color;
axis equal,  xlabel('x [AU]'), ylabel('y [AU]'), zlabel('z [AU]')


end