function [THS, alphaHS, betaHS, XHS, XSOL,fval,exitflag,output,lambda,grad,hessian] = runDirectTranscription(t0, TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, href, hh, v1, v2, muS, data, Bounds)
%% DIRECT TRANSCRIPTION

[ m, T, r, z, ~, vr, ~, vz, ~, ~, ~, TH, ~, ~, ~, gamma, ~, ~, ~, ~, ~, ~, T_inplane, T_outplane, theta_dot, time, ~] = ...
    Conway(TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, href, hh, v1, v2, muS, data);

X = [r', TH', z', vr', theta_dot', vz', m'];
options = optimset( 'Algorithm','interior-point', ...
                        'Display','iter', ... %          
                        'LargeScale','on', ...
                        'MaxIter',50000, ...
                        'MaxFunEvals', 5000000);

N = data.n_int;
DU = astroConstants(2);
TU = (DU^3/muS).^0.5;
MU = data.Mdry;


XX0 = zeros(1, N*10);
Xad = zeros(N,7); Xd = Xad;
%adimensionalization of initial conditions
Xad(:,1) = X(:,1)/DU;
Xad(:,2) = X(:,2); %already adimensional
Xad(:,3) = X(:,3)/DU;
Xad(:,4) = X(:,4)/DU*TU; 
Xad(:,5) = X(:,5)*TU; 
Xad(:,6) = X(:,6)/DU*TU;
Xad(:,7) = X(:,7)/MU;

% thr = 1e-3;
% ubeta = atan2((T_outplane < thr).*thr + (T_outplane > thr).*T_outplane, T_inplane);
% ubeta = atan2(T_outplane, T_inplane);
ubeta = asin(T_outplane./ T);
ualpha = gamma;

data.xi = Xad(1,:);
data.xf = Xad(end,:);
data.muS = muS;   
timead = time/TU;
data.time = timead;

%%%%%

Xpropad = zeros(N, 7);
Xpropad(1,:) = Xad(1,:); %first value equal

  %RK4 forward integration (in each interval of delta t)
    for k = 1:N-1
        xx = Xpropad(k,:);

        hhh = timead(k+1) - timead(k);
        K1 = EoMpolarAD(timead(k), xx, T(k), ualpha(k), ubeta(k), muS, data);
        K2 = EoMpolarAD(timead(k)+ hhh/2, xx + hhh/2*K1, 0.5*(T(k+1) + T(k)), 0.5*(ualpha(k+1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)), muS, data);
        K3 = EoMpolarAD(timead(k)+ hhh/2, xx + hhh/2*K2, 0.5*(T(k+1) + T(k)), 0.5*(ualpha(k+1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)),  muS, data);
        K4 = EoMpolarAD(timead(k)+ hhh, xx + hhh*K3, T(k+1), ualpha(k+1), ubeta(k+1), muS, data);    

        Xpropad(k+1,:) = Xpropad(k,:) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end
Xd(:,1) = Xpropad(:,1)*DU;
Xd(:,2) = Xpropad(:,2); %already adimensional
Xd(:,3) = Xpropad(:,3)*DU;
Xd(:,4) = Xpropad(:,4)*DU/TU; 
Xd(:,5) = Xpropad(:,5)/TU; 
Xd(:,6) = Xpropad(:,6)*DU/TU;
Xd(:,7) = Xpropad(:,7)*MU;

figure()
sgtitle('confrontatio of propagated and adimensional states..')
subplot(7,1,1), plot(Xad(:,1)), hold on, plot(Xpropad(:,1))
subplot(7,1,2), plot(Xad(:,2)), hold on, plot(Xpropad(:,2))
subplot(7,1,3), plot(Xad(:,3)), hold on, plot(Xpropad(:,3))
subplot(7,1,4), plot(Xad(:,4)), hold on, plot(Xpropad(:,4))
subplot(7,1,5), plot(Xad(:,5)), hold on, plot(Xpropad(:,5))
subplot(7,1,6), plot(Xad(:,6)), hold on, plot(Xpropad(:,6))
subplot(7,1,7), plot(Xad(:,7)), hold on, plot(Xpropad(:,7)), legend('ad','prop')
%%%%%
fprintf('Press any button to continue...')
pause()

%already adimensional
data.xi = Xad(1,:);
data.xf = Xad(end,:);

% Adimensionalization of variables
figure(),
sgtitle('Conway Control'),
subplot(4,2,[1 2]), plot(timead, rad2deg(ubeta),'HandleVisibility','off'), hold on, xlabel('time [TU]')
yline(rad2deg(Bounds.beta_ub),'r','DisplayName','UB'); hold on, yline(rad2deg(Bounds.beta_lb),'b','DisplayName','LB'),
title('$\beta$'), grid on, legend()
subplot(4,2,[3 4]), plot(timead, rad2deg(gamma),'HandleVisibility','off'), hold on, xlabel('time [TU]')
yline(rad2deg(Bounds.alpha_ub),'r','DisplayName','UB'); hold on, yline(rad2deg(Bounds.alpha_lb),'b','DisplayName','LB'),
title('$\alpha$'), grid on, legend()
subplot(4,2,[5 6]),  plot(timead, T,'HandleVisibility','off'), title('T'), grid on, xlabel('time [TU]')
yline(data.Tmax,'r','DisplayName','UB'); hold on, yline(0,'b','DisplayName','LB'), legend()
subplot(4,2,7), plot(timead, T_inplane,'HandleVisibility','off'), title('T inplane'), grid on, xlabel('time [TU]')
subplot(4,2,8),  plot(timead, T_outplane,'HandleVisibility','off'), title('T outplane'), grid on, xlabel('time [TU]')
 
%definition of the initial guess (sub-optimal conway solution)
for ii = 1:N
    XX0((ii-1)*10 +1: (ii-1)*10+10) = [Xad(ii,:) T(ii) ualpha(ii) ubeta(ii)];
end
       
fprintf('check on the aiming of mars and boundaries')
kepMtry = uplanet(t0 + TOF,4);
rMtry = kep2car2(kepMtry, muS);
R = refplane2car( X(end,1), X(end,3),  X(end,1)*X(end,5), X(end,4), X(end,6), X(end,2), r1vers, href);
rMtry - R
%definition of upper boundary and lower boundary
[LB, UB] = LBUB(XX0, Xad,  data, Bounds);
find(XX0 < LB)
find(XX0 > UB)
fprintf('Press any button to continue...')
pause()
[XSOL,fval,exitflag,output,lambda,grad,hessian] = fmincon(@DTmethod, XX0,[],[],[],[],LB,UB,@EoMconstr, options, data);
%%
XHS = zeros(N, 7);
THS = zeros(1,N);
alphaHS = THS;
betaHS = THS;
%reconstruction of states
for ii = 1:length(TH)
    [XHS(ii,:)] = XSOL((ii-1)*10 +1: (ii-1)*10+7);
     XHS(ii,1) = XHS(ii,1) * DU;
     XHS(ii,2) = XHS(ii,2);
     XHS(ii,3) = XHS(ii,3) * DU;
     XHS(ii,4) = XHS(ii,4) * DU/TU;
     XHS(ii,5) = XHS(ii,5) / TU;
     XHS(ii,6) = XHS(ii,6) * DU/TU;
     XHS(ii,7) = XHS(ii,7) * MU;
    [THS(ii)] = XSOL((ii-1)*10+8);
    [alphaHS(ii)] = XSOL((ii-1)*10+9);
    [betaHS(ii)] = XSOL((ii-1)*10+10);
end
 
rHS = XHS(:,1); thHS = XHS(:,2); zHS = XHS(:,3); vrHS = XHS(:,4);
thdHS = XHS(:,5); vzHS = XHS(:,6); mHS = XHS(:,7); 

sHS = (rHS.^2 + zHS.^2).^0.5;
maxT = d2T (sHS/DU);
maxT = 0.25*(maxT >= 0.25) + maxT.*(maxT < 0.25);


figure()
sgtitle('Optimal Control')
subplot(3,1,1), plot(timead, THS, 'DisplayName', 'Optimal control'), hold on, plot(timead, T,  'DisplayName', 'Conway Solution'), hold on, 
plot(timead, maxT, 'b:','DisplayName','max available T'), title('T'), ylabel('T [N]'), xlabel('t [TU]'), legend()
subplot(3,1,2), plot(timead, rad2deg(alphaHS),'DisplayName', 'Optimal control'), hold on, plot(timead, rad2deg(ualpha),  'DisplayName', 'Conway Solution'), hold on, title('$\alpha$'), ylabel('$\alpha$ [deg]'), xlabel('t [TU]'), legend()
subplot(3,1,3), plot(timead, rad2deg(betaHS),'DisplayName', 'Optimal control'), hold on, plot(timead, rad2deg(ubeta),  'DisplayName', 'Conway Solution'), hold on, title('$\beta$'), ylabel('$\beta$ [deg]'), xlabel('t [TU]'), legend()


figure()
subplot(4,2,1), plot(timead, rHS,'DisplayName', 'Optimal control'), hold on, plot(timead, r,  'DisplayName', 'Conway Solution'), hold on, ylabel('r [km]')
subplot(4,2,2), plot(timead, thHS,'DisplayName', 'Optimal control'), hold on, plot(timead, TH,  'DisplayName', 'Conway Solution'), hold on,legend(), ylabel('$\theta$ [rad]')
subplot(4,2,3), plot(timead, zHS,'DisplayName', 'Optimal control'), hold on, plot(timead, z,  'DisplayName', 'Conway Solution'), hold on,ylabel('z [km]')
subplot(4,2,4), plot(timead, vrHS,'DisplayName', 'Optimal control'), hold on, plot(timead, vr,  'DisplayName', 'Conway Solution'), hold on,ylabel('$v_r$ [km/s]')
subplot(4,2,5), plot(timead, thdHS,'DisplayName', 'Optimal control'), hold on, plot(timead, theta_dot,  'DisplayName', 'Conway Solution'), hold on,ylabel('$\theta_dot$ [km]')
subplot(4,2,6), plot(timead, vzHS,'DisplayName', 'Optimal control'), hold on, plot(timead, vz,  'DisplayName', 'Conway Solution'), hold on,ylabel('$v_z$ [km/s]')
subplot(4,2,7), plot(timead, mHS,'DisplayName', 'Optimal control'), hold on, plot(timead, m,  'DisplayName', 'Conway Solution'), hold on,ylabel('m [kg]')
subplot(4,2,8), plot(timead, THS,'DisplayName', 'Optimal control'), hold on, plot(timead, T,  'HandleVisibility','off'), hold on,ylabel('T [N]')
plot(timead, maxT, 'b:','DisplayName','max available T')

end

