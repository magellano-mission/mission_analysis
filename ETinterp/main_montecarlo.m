%% MONTECARLO SIMULATION
%giving the same control law (no predictive control) what happens to the
%state...
% clearvars -except dataNS1 dataNS2 dataNS3 dataRS1 dataRS2 dataECS

%%
close all, clc
data = dataNS1;

NMC = 100;
X0 = [r1(1); TH1(1); z1(1); vr1(1); theta_dot1(1); vz1(1); m1(1)]';
XMC = zeros(data.n_int, 7,NMC); 

alpha = gamma1;
beta = asin(data.comps.T_out./ T1);

betaMC = zeros(data.n_int,NMC);
alphaMC = betaMC; TMC = betaMC;

sigma_T   = 0*max(abs(T1));
sigma_alpha   = deg2rad(1); %1e-1*max(abs(alpha));
sigma_beta    = deg2rad(1); %1e-1*max(abs(beta));

fprintf('sigma_T \t %d deg \n', sigma_T );
fprintf('sigma_alpha \t %d deg \n', rad2deg(sigma_alpha) );
fprintf('sigma_beta \t %d deg \n', rad2deg(sigma_beta));

sigma_inj = 0*1e-1;

wb = waitbar(0,'MONTECARLO SIMULATION');
vr = zeros(data.n_int, NMC);
vt = vr; r = vr; th = r; z = r; vz = vt;

for nmc = 1:NMC
   waitbar(nmc/NMC, wb);
   err = sigma_inj*randn(1,7); err(7) = 0;
   X0MC = X0 + 0*err;
   
   %introducing accuracy
   alphaMC(:,nmc)    = alpha + sigma_alpha*randn(1,data.n_int);
   betaMC(:,nmc)     = beta + sigma_beta*randn(1,data.n_int);
   TMC(:,nmc)        = T1 + sigma_T*randn(1,data.n_int);
   
   %propagation
  [XMC(:,:,nmc)] = EoMpropRK4_2(X0MC, time1, TMC(:,nmc), alphaMC(:,nmc), betaMC(:,nmc), muS, data, 0);
  r(:,nmc) = squeeze(XMC(:,1,nmc));
  th(:,nmc)  = squeeze(XMC(:,2,nmc));
  z(:,nmc)   = squeeze(XMC(:,3,nmc));
  vr(:,nmc)  = squeeze(XMC(:,4,nmc));
  vt(:,nmc)  = r(:,nmc).*squeeze(XMC(:,5,nmc));
  vz(:,nmc)  = squeeze(XMC(:,6,nmc));
end   
delete(wb)

rrMC = zeros(data.n_int, 3, NMC);vvMC = rrMC;
%computation of the state of 
for nmc = 1:NMC
    for i = 1:data.n_int
        [rrMC(i,:,nmc), vvMC(i,:,nmc)] = refplane2car( r(i,nmc), z(i,nmc),  vt(i,nmc), vr(i,nmc), vz(i,nmc), th(i,nmc), r1vers_1, hvers_1);
    end
end

  
%stats
%statistics 
accuracy = zeros(data.n_int, 3);
mean = zeros(data.n_int,3);

for t = 1:data.n_int
mean(t,1) = sum(alphaMC(t,:))/NMC;
mean(t,2) = sum(betaMC(t,:))/NMC;
mean(t,3) = sum(TMC(t,:))/NMC;

accuracy(t,1) = sqrt(sum((alphaMC(t,:) - mean(t,1)).^2)/(NMC-1));
accuracy(t,2) = sqrt(sum((betaMC(t,:) - mean(t,2)).^2)/(NMC-1));
accuracy(t,3) = sqrt(sum((TMC(t,:) - mean(t,3)).^2)/(NMC-1));

end

figure()
RE = zeros(length(TOFr1), 3); RM = RE;
for i =1:length(TOFr1)
    kepE = uplanet(t01+TOFr1(i), 3);
    kepM = uplanet(t01+TOFr1(i), 4);
    RE(i,:) = kep2car2(kepE, muS);
    RM(i,:) = kep2car2(kepM, muS);
end
% clear kepE kepM
mpp = plot3(RM(:,1),RM(:,2),RM(:,3),'HandleVisibility','Off'); hold on
epp = plot3(RE(:,1),RE(:,2),RE(:,3),'HandleVisibility','Off'); hold on
plot3(rrNS1(:,1), rrNS1(:,2), rrNS1(:,3),'m', 'DisplayName','Conway Solution'), hold on
sgtitle('Injection maneuver error')
for nmc = 1:NMC
    rrMCtemp = squeeze(rrMC(:,:,nmc));
    ppp = plot3(rrMCtemp(:,1),rrMCtemp(:,2),rrMCtemp(:,3),'k:','HandleVisibility','Off'); hold on
    ppp.Color(4) = 0.2;
end

ppp.HandleVisibility = 'on';
ppp.DisplayName = 'Montecarlo Simulation';

mpp0 = plot3(RM(end,1),RM(end,2),RM(end,3), 'o','MarkerSize',10,'DisplayName','Mars @ Arrival'); hold on
epp0 = plot3(RE(1,1),RE(1,2),RE(1,3),'o','MarkerSize',10,'DisplayName','Earth @ Departure'); hold on
mpp0.Color = mpp.Color;
epp0.Color = epp.Color;
plotSun(), legend()

figure()
for kk = 1:NMC
    
    subplot(2,3,1), 
    gammaplot = plot(TOFr1, rad2deg(alphaMC(:,kk)),'k','HandleVisibility','off'); hold on
    gammaplot.Color(4) = 0.3;
    
    subplot(2,3,2), 
    betaplot = plot(TOFr1, rad2deg(betaMC(:,kk)),'k','HandleVisibility','off'); hold on
    betaplot.Color(4) = 0.3;
    
    subplot(2,3,3), 
    Tplot = plot(TOFr1, TMC(:,kk),'k','HandleVisibility','off'); hold on
    Tplot.Color(4) = 0.3;

end

Tplot.HandleVisibility = 'on';
Tplot.DisplayName = 'T Montecarlo Simulation';

betaplot.HandleVisibility = 'on';
betaplot.DisplayName = '$\beta$ Montecarlo Simulation';

gammaplot.HandleVisibility = 'on';
gammaplot.DisplayName = '$\alpha$ Montecarlo Simulation';

subplot(2,3,1), plot(TOFr1, rad2deg(mean(:,1)), 'DisplayName','Mean Value'), 
ylabel('deg'), xlabel('T days')
title('$\alpha$'), legend()
subplot(2,3,2), title('$\alpha$'), legend()

subplot(2,3,2), plot(TOFr1, rad2deg(mean(:,2)), 'DisplayName','Mean Value'), 
ylabel('deg'), xlabel('T days')
title('$\beta$'), legend()
subplot(2,3,2), title('$\beta$'), legend()

subplot(2,3,3), plot(TOFr1, mean(:,3), 'DisplayName','Mean Value'), 
ylabel('N'), xlabel('T days')
title('$T$'), legend()

subplot(2,3,4),
plot(TOFr1, rad2deg(accuracy(:,1))), ylabel('deg'), xlabel('T days')
title('accuracy $\sigma_{\alpha}$')

subplot(2,3,5),
plot(TOFr1, rad2deg(accuracy(:,2))), ylabel('deg'), xlabel('T days')
title('accuracy $\sigma_{\beta}$')

subplot(2,3,6),
plot(TOFr1, accuracy(:,3)), ylabel('N'), xlabel('T days')
title('accuracy $\sigma_{T}$')
