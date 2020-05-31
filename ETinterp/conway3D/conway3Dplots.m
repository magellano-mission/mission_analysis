%conway 3d plots
DU = astroConstants(2);
TU = (DU^3/muS)^0.5;

if ~isnan(r)   
    
[kepEarth, muS] = uplanet(t0, 3);
[rE0, vE0] = kep2car2(kepEarth, muS);
[kepMars, muS] = uplanet(t0 + TOF, 4);
[rMend, vMend] = kep2car2(kepMars, muS);
  
r1vers = rE0/norm(rE0);
r2vers = rMend/norm(rMend) ; 
    
R1 = zeros(length(TOFr), 3); RM = R1;
REnorm = zeros(length(TOFr),1); RMnorm = REnorm;
rr = R1; vv = rr;
for i =1:length(TOFr)
    kepE = uplanet(t0+TOFr(i), 3);
    kepM = uplanet(t0+TOFr(i), 4);
    R1(i,:) = kep2car2(kepE, muS);
    RM(i,:) = kep2car2(kepM, muS);
    REnorm(i) = norm(R1(i,:));
    RMnorm(i) = norm(RM(i,:));
    [rr(i,:), vv(i,:)] = refplane2car( r(i), z(i),  vt(i), vr(i), TH(i), r1vers, href);
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
fprintf('Isp: \t %d s\n', data_stacks.Isp)
fprintf('Mdry: \t %d kg\n', data_stacks.Mdry)
fprintf('Propellant mass: \t %d kg \n', m(1) - m(end))
fprintf('Mass ratio: \t %d \n', m(end)/m(1))
fprintf('Fuel Mass Fraction: \t %d \n', (m(1) - m(end))/m(1))
fprintf('max T : \t %d N\n',max(abs(T)))

figure()
sgtitle('')
% subplot(2,2,[1 2])
rmp = plot3(RM(:,1)/DU, RM(:,2)/DU, RM(:,3)/DU,'HandleVisibility','off'); hold on,
rep = plot3(R1(:,1)/DU, R1(:,2)/DU, R1(:,3)/DU,'HandleVisibility','off'); hold on,
rrp = plot3(rr(:,1)/DU, rr(:,2)/DU, rr(:,3)/DU,'--','HandleVisibility','off'); hold on,
rmpend = plot3(RM(end,1)/DU, RM(end,2)/DU, RM(end,3)/DU,'o','DisplayName','Mars @ Arrival'); hold on,
rep1 = plot3(R1(1,1)/DU, R1(1,2)/DU, R1(1,3)/DU,'o','DisplayName','Earth @ Departure'); hold on,
rrp1 = plot3(rr(1,1)/DU, rr(1,2)/DU, rr(1,3)/DU,'+','HandleVisibility','off'); hold on,
rrpend = plot3(rr(end,1)/DU, rr(end,2)/DU, rr(end,3)/DU,'+','HandleVisibility','off');hold on
I = imread('Sun.jpg'); RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];  RI.YWorldLimits = [-90 90]; 
rSun = 20*almanac('Sun','Radius','kilometers','sphere');
[XSun, YSun, ZSun] = ellipsoid(0, 0, 0, rSun, rSun, rSun, 100); % spheric centered Mars
planet = surf(XSun/DU, YSun/DU, -ZSun/DU,'Edgecolor', 'none','HandleVisibility','off'); hold on
set(planet,'FaceColor','texturemap','Cdata',I), axis equal
legend()
rmpend.Color = rmp.Color; rep1.Color = rep.Color; rrp1.Color = rrp.Color; rrpend.Color = rrp.Color;
axis equal,  xlabel('x [AU]'), ylabel('y [AU]'), zlabel('z [AU]')

figure()
subplot(2,1,1), 
cM = plot(TOFr, RMnorm/DU,'Displayname','$r_{Mars}$'); hold on, 
cE = plot(TOFr, REnorm/DU, 'Displayname','$r_{Earth}$'); hold on, 
cS = plot(TOFr, r/DU,'Displayname','$r_{Mars}$');  hold off, 
title('In-plane motion'), xlabel('TOF [days]'), ylabel('r [AU]')
subplot(2,1,2), 
plot(TOFr, RM*href/DU), hold on, plot(TOFr, R1*href/DU), hold on, plot(TOFr, z/DU), 
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
vr = plot(TOFr, vr, 'DisplayName','$v_{r}$'); hold on,
vte.Color = cE.Color;
vtm.Color = cM.Color;
vsn.Color = cS.Color;
vts.Color = cS.Color;
legend(), title('Velocity Components')


else
    fprintf('No real solution for Conway algorithm \n')
end
 