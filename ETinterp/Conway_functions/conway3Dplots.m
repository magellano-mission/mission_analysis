function [rr, vv, T_cart] = conway3Dplots(t0, TOF, href, N_rev, q, m, T, r, z, vr, vt, vz, acc_inplane, acc_out, acc, TH, gamma1, gamma2, gamma, v1perp, v2perp, v1tra, v2tra, vnorm, T_inplane, T_outplane, TOFr, data)
%output: cartesian components of states
%conway 3d plots
if ~isnan(r)   
    
[kepEarth, muS] = uplanet(t0, 3);
[rE0, vE0] = kep2car2(kepEarth, muS);
[kepMars, muS] = uplanet(t0 + TOF, 4);
[rMend, vMend] = kep2car2(kepMars, muS);

AU = astroConstants(2);
TU = (AU^3/muS).^0.5;

r1vers = rE0/norm(rE0);
r2vers = rMend/norm(rMend) ; 
    
R1 = zeros(length(TOFr), 3); RM = R1;
REnorm = zeros(length(TOFr),1); RMnorm = REnorm;
rr = R1; vv = rr; T_cart = rr;
for i =1:length(TOFr)
    kepE = uplanet(t0+TOFr(i), 3);
    kepM = uplanet(t0+TOFr(i), 4);
    R1(i,:) = kep2car2(kepE, muS);
    RM(i,:) = kep2car2(kepM, muS);
    REnorm(i) = norm(R1(i,:));
    RMnorm(i) = norm(RM(i,:));
    [rr(i,:), vv(i,:)] = refplane2car(r(i), z(i),  vt(i), vr(i), vz(i), TH(i), r1vers, href);
    [T_cart(i,:)] = thrust_cart( T_inplane(i), T_outplane(i), gamma(i), TH(i), r1vers, href);
end

dirt1 = cross(href , r1vers);
dirt2 = cross(href , r2vers);
 
vr1 = r1vers * vr(1);
vr2 = r2vers * vr(end);
 
vt1 = dirt1 * vt(1);
vt2 = dirt2 * vt(end);
 
v1 = vr1 + vt1 + v1perp*href; v2 = vr2 + vt2 + v2perp*href;
v_inf1 = v1 - vE0; v_inf2 = vMend - v2;
 
% optimization result
datedep=mjd20002date(t0);
datearr=mjd20002date(t0 + TOF);
fprintf('OPTIMIZATION RESULTS \n')
fprintf('Departure Date \t %d %d %d \n', datedep(3), datedep(2),datedep(1));
fprintf('Arrival Date \t %d %d %d \n', datearr(3), datearr(2),datearr(1));
fprintf('TOF \t \t %d days (%d yrs)\n', TOF, TOF/365);
fprintf('N_rev \t \t %d \n', N_rev);
fprintf('q \t \t %d \n \n', q);
fprintf('Departure v_inf: \t %d km/s (C3 = %d km^2/s^2) \n', norm(v_inf1), norm(v_inf1)^2)
fprintf('Arrival v_inf: \t %d km/s (C3 = %d km^2/s^2)   \n \n', norm(v_inf2), norm(v_inf2)^2)
% fprintf('Isp: \t %d s\n', data.Isp)
fprintf('Mdry: \t %d kg\n', data.Mdry)
fprintf('Propellant mass: \t %d kg \n', m(1) - m(end))
fprintf('Mass ratio: \t %d \n', m(end)/m(1))
fprintf('Fuel Mass Fraction: \t %d \n', (m(1) - m(end))/m(1))
fprintf('max T : \t %d N\n',max(abs(T)))

figure()
sgtitle('')
% subplot(2,2,[1 2])
rmp = plot3(RM(:,1), RM(:,2), RM(:,3),'HandleVisibility','off'); hold on,
rep = plot3(R1(:,1), R1(:,2), R1(:,3),'HandleVisibility','off'); hold on,
rrp = plot3(rr(:,1), rr(:,2), rr(:,3),'--','HandleVisibility','off'); hold on,
rmpend = plot3(RM(end,1), RM(end,2), RM(end,3),'o','DisplayName','Mars @ Arrival'); hold on,
rep1 = plot3(R1(1,1), R1(1,2), R1(1,3),'o','DisplayName','Earth @ Departure'); hold on,
rrp1 = plot3(rr(1,1), rr(1,2), rr(1,3),'+','HandleVisibility','off'); hold on,
rrpend = plot3(rr(end,1), rr(end,2), rr(end,3),'+','HandleVisibility','off');hold on

%velocity plot
for i = 1:20:length(TOFr)
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


figure()
subplot(2,1,1), 
cM = plot(TOFr, RMnorm/AU,'Displayname','$r_{Mars}$'); hold on, 
cE = plot(TOFr, REnorm/AU, 'Displayname','$r_{Earth}$'); hold on, 
cS = plot(TOFr, r/AU,'Displayname','$r_{Mars}$');  hold off, 
title('In-plane motion'), xlabel('TOF [days]'), ylabel('r [AU]')
subplot(2,1,2), 
plot(TOFr, RM*href/AU), hold on, plot(TOFr, R1*href/AU), hold on, plot(TOFr, z/AU), 
hold off, title('Out-of-plane motion'), xlabel('TOF [days]'), ylabel('z [AU]')

figure()
sgtitle('Thrust Profile')
subplot(5,2,1), plot(TOFr, acc_inplane), title('$a_{inplane}$')
subplot(5,2,2), plot(TOFr, acc_out), title('$a_{outplane}$')
subplot(5,2,3), plot(TOFr, T_inplane), title('$T_{inplane}$')
subplot(5,2,4), plot(TOFr, T_outplane), title('$T_{outplane}$')
subplot(5,2,[5 6]), plot(TOFr, acc), title('$a_{tot}$')
subplot(5,2,[7 8]), plot(TOFr, T), title('$T$')
subplot(5,2,[9 10]), plot(TOFr, m), title('$m$')

figure()
sgtitle('Boundary Conditions')
subplot(2,1,1),
gs = plot(TOFr, rad2deg(gamma),'DisplayName','spacecraft'); hold on 
ge = yline(rad2deg(gamma1),'DisplayName','Earth Departure'); hold on, 
gm = yline(rad2deg(gamma2),'DisplayName','Mars Arrival'); hold off
gs.Color = cS.Color;
gm.Color = cM.Color;
ge.Color = cE.Color;
legend(), title('Flight path angle')


subplot(2,1,2)
vts = plot(TOFr, vt,'DisplayName','$v_{t}$'); hold on
vte = yline(norm(v1tra),'DisplayName','$v^E_\theta$'); hold on, 
vtm = yline(norm(v2tra),'DisplayName','$v^M_\theta$'); hold on,
vsn = plot(TOFr, vnorm,'DisplayName','$|v|$'); hold on, 
vrplot = plot(TOFr, vr, 'DisplayName','$v_{r}$'); hold on,
vte.Color = cE.Color;
vtm.Color = cM.Color;
vsn.Color = cS.Color;
vts.Color = cS.Color;
legend(), title('Velocity Components')


else
    fprintf('No real solution for Conway algorithm \n')
end
end