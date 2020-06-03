%% MONTECARLO SIMULATION
%giving the same control law (no predictive control) what happens to the
%state...
% MISSINIG INITALIZATION....
NMC = 100;
X0 = [r(1), TH(1), z(1), vr(1), theta_dot(1), vz(1), m(1)];
XMC = zeros(data.n_int, 7,NMC); 
sigma_inj = X0*0.001.*ones(1,7); 
sigma_inj(7) = 0; %same mass of propellant of course...

wb = waitbar(0,'MONTECARLO SIMULATION');
for nmc = 1:NMC
   waitbar(nmc/NMC, wb);
   X0MC = X0 + sigma_inj.*randn(1,7);
  [XMC(:,:,nmc)] = EoMpropRK4(X0MC, time1, acc_inplane1, acc_out1, muS, data, 0);
end   
delete wb
rrMC = zeros(data.n_int, 3, NMC);vvMC = rrMC;
pause()
%computation of the state of 
for nmc = 1:NMC
    r   = squeeze(XMC(:,1,nmc));
    th  = squeeze(XMC(:,2,nmc));
    z   = squeeze(XMC(:,3,nmc));
    vr  = squeeze(XMC(:,4,nmc));
    vt  = r.*squeeze(XMC(:,5,nmc));
    vz  = squeeze(XMC(:,6,nmc));
    for i = 1:data.n_int
        [rrMC(i,:,nmc), vvMC(i,:,nmc)] = refplane2car( r(i), z(i),  vt(i), vr(i), vz(i), th(i), r1vers_1, hvers_1);
    end
end
%%
figure()
RE = zeros(length(TOFr1), 3); RM = RE;
for i =1:length(TOFr1)
    kepE = uplanet(t01+TOFr1(i), 3);
    kepM = uplanet(t01+TOFr1(i), 4);
    RE(i,:) = kep2car2(kepE, muS);
    RM(i,:) = kep2car2(kepM, muS);
end
clear kepE kepM
mpp = plot3(RM(:,1),RM(:,2),RM(:,3),'HandleVisibility','Off'); hold on
epp = plot3(RE(:,1),RE(:,2),RE(:,3),'HandleVisibility','Off'); hold on
plot3(rrNS1(:,1), rrNS1(:,2), rrNS1(:,3),'m', 'LineWidth',5, 'DisplayName','Conway Solution'), hold on
sgtitle('Injection maneuver error')
for nmc = 1:NMC
    ppp = plot3(rrMC(:,1,nmc),rrMC(:,2,nmc),rrMC(:,3,nmc),'k:','HandleVisibility','Off'); hold on
    ppp.Color(4) = 0.2;
end
ppp.HandleVisibility = 'on';
ppp.DisplayName = 'Montecarlo Simulation';

mpp0 = plot3(RM(end,1),RM(end,2),RM(end,3), 'o','MarkerSize',10,'DisplayName','Mars @ Arrival'); hold on
epp0 = plot3(RE(1,1),RE(1,2),RE(1,3),'o','MarkerSize',10,'DisplayName','Earth @ Departure'); hold on
mpp0.Color = mpp.Color;
epp0.Color = epp.Color;
plotSun(), legend()