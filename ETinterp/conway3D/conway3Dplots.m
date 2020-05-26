%conway 3d plots

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
    [rr(i,:), vv(i,:)] = refplane2car( r(i), z(i),  vt(i), vr(i), TH(i), r1vers, RCRRv);
end
 
 
dirt1 = cross(RCRRv , r1vers);
dirt2 = cross(RCRRv , r2vers);
 
vr1 = r1vers * vr(1);
vr2 = r2vers * vr(end);
 
vt1 = dirt1 * vt(1);
vt2 = dirt2 * vt(end);
 
v1 = vr1 + vt1 + v1perp*RCRRv; v2 = vr2 + vt2 + v2perp*RCRRv;
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
subplot(2,2,[1 2])
rmp = plot3(RM(:,1), RM(:,2), RM(:,3),'HandleVisibility','off'); hold on,
rep = plot3(R1(:,1), R1(:,2), R1(:,3),'HandleVisibility','off'); hold on,
rrp = plot3(rr(:,1), rr(:,2), rr(:,3),'--','HandleVisibility','off'); hold on,
rmpend = plot3(RM(end,1), RM(end,2), RM(end,3),'o','DisplayName','Mars @ Arrival'); hold on,
rep1 = plot3(R1(1,1), R1(1,2), R1(1,3),'o','DisplayName','Earth @ Departure'); hold on,
rrp1 = plot3(rr(1,1), rr(1,2), rr(1,3),'+','HandleVisibility','off'); hold on,
rrpend = plot3(rr(end,1), rr(end,2), rr(end,3),'+','HandleVisibility','off'); hold on, legend()

rmpend.Color = rmp.Color; rep1.Color = rep.Color; rrp1.Color = rrp.Color; rrpend.Color = rrp.Color;
axis equal,  title('complete path')
 
subplot(2,2,3), 
plot(TOFr, RMnorm), hold on, plot(TOFr, REnorm), hold on, plot(TOFr, r), 
hold off, title('in-plane motion')
subplot(2,2,4), 
plot(TOFr, RM*RCRRv), hold on, plot(TOFr, R1*RCRRv), hold on, plot(TOFr, z), 
hold off, title('out-of-plane motion')

figure()
sgtitle('Thrust Profile')
subplot(5,2,1), plot(TOFr, acc_inplane), title('a_{inplane}')
subplot(5,2,2), plot(TOFr, acc_out), title('a_{outplane}')
subplot(5,2,3), plot(TOFr, T_inplane), title('T_{inplane}')
subplot(5,2,4), plot(TOFr, T_outplane), title('T_{outplane}')
subplot(5,2,[5 6]), plot(TOFr, acc), title('a_{tot}')
subplot(5,2,[7 8]), plot(TOFr, T), title('T')
subplot(5,2,[9 10]), plot(TOFr, m), title('m')

figure()
sgtitle('Flight path angle')
plot(TOFr, gamma,'DisplayName','spacecraft'), hold on 
yline(gamma1,'DisplayName','Departure'); hold on, yline(gamma2,'DisplayName','Mars');
legend(),

figure()
plot(TOFr, vt,'DisplayName','$v_{t}$'), hold on
yline(norm(v1tra),'DisplayName','$v^E_\theta$'); hold on, yline(norm(v2tra),'DisplayName','$v^M_\theta$');
plot(TOFr, vnorm,'DisplayName','$|v|$'), hold on, plot(TOFr, vr, 'DisplayName','$v_{r}$')
legend()
%  annotation(gcf,'textbox',...
%     [0.15 0.75 0.29 0.15],...
%     'String',{strcat('Dep :',num2str(datedep(3)),'.', num2str(datedep(2)),'.', num2str(datedep(1)), '...', ...
%                      'Arr :', num2str(datearr(3)),'.', num2str(datearr(2)),'.', num2str(datearr(1))), ...
%               strcat('Nrev :',num2str(N_rev),'...','TOF :',num2str(TOF),'days'), ...
%               strcat('v_{inf} dep :', num2str(v_inf1),'km/s', 'v_{inf} arr :',num2str(v_inf2), 'km/s'), ...
%               strcat('max|T| :', num2str(max(abs(T))),'N','...','Isp :',num2str(data_stacks.Isp),'s','...','Mdry :', num2str(data_stacks.Mdry),'kg')},...
%     'FitBoxToText','on');
else
    fprintf('No real solution for Conway algorithm \n')
end
 