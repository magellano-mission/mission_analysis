%direct transcription

options = optimoptions('fmincon','Display','iter', ... %                        'MaxFunctionEvaluations',1e4, ...
                        'Algorithm','interior-point', ...
                        'MaxIter',50000, ...
                        'MaxFunEvals', 500000);

data = data_stacks;
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

ubeta = atan(T_outplane./T_inplane);
ualpha = gamma;

data.xi = Xad(1,:);
data.xf = Xad(end,:);
data.muS = muS;   

Xpropad = zeros(N, 7);
Xpropad(1,:) = Xad(1,:); %first value equal
timead = time/TU;
data.time = timead;
    %RK4 forward integration (in each interval of delta t)
    for k = 1:N-1
        xx = Xpropad(k,:);
%         xx = x(k,:);

        hhh = timead(k+1) - timead(k);
        K1 = polar_ad(timead(k), xx, T(k), ualpha(k), ubeta(k), muS, data);
        K2 = polar_ad(timead(k)+ hhh/2, xx + hhh/2*K1, 0.5*(T(k+1) + T(k)), 0.5*(ualpha(k+1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)), muS, data);
        K3 = polar_ad(timead(k)+ hhh/2, xx + hhh/2*K2, 0.5*(T(k+1) + T(k)), 0.5*(ualpha(k+1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)),  muS, data);
        K4 = polar_ad(timead(k)+ hhh, xx + hhh*K3, T(k+1), ualpha(k+1), ubeta(k+1), muS, data);    

%         Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
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
subplot(7,2,1), plot(Xad(:,1)-Xpropad(:,1))
subplot(7,2,3), plot(Xad(:,2)-Xpropad(:,2))
subplot(7,2,5), plot(Xad(:,3)-Xpropad(:,3))
subplot(7,2,7), plot(Xad(:,4)-Xpropad(:,4))
subplot(7,2,9), plot(Xad(:,5)-Xpropad(:,5))
subplot(7,2,11), plot(Xad(:,6)-Xpropad(:,6))
subplot(7,2,13), plot(Xad(:,7)-Xpropad(:,7))

%already adimensional
data.xi = Xpropad(1,:);
data.xf = Xpropad(end,:);
%% Adimensionalization of variables
ubeta = atan( max( min( (T_outplane )./(T_inplane), 1) , -1) );

figure(),subplot(5,1,1), plot(rad2deg(ubeta)), title('Beta'), grid on
subplot(5,1,2), plot(T_inplane), title('T inplane'), grid on
subplot(5,1,3), plot(T_outplane), title('T outplane'), grid on
subplot(5,1,4), plot(T), title('T'), grid on
ualpha = gamma;
%definition of the initial guess (sub-optimal conway solution)
for ii = 1:N
    XX0((ii-1)*10 +1: (ii-1)*10+10) = [Xpropad(ii,:) T(ii) ualpha(ii) ubeta(ii)];
end
%%
%definition of upper boundary and lower boundary
[LB, UB] = LBUB(XX0, data);

find(XX0 < LB)
% [c, ceq] = EoM(XX0, data);
%%
[XSOL,fval,exitflag,output,lambda,grad,hessian] = fmincon(@DTmethod, XX0,[],[],[],[],LB,UB,@EoM, options, data);
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

rHS = XHS(:,1);
figure()
sgtitle('Control allocation')
subplot(3,1,1), plot(timead, THS), hold on, plot(timead, T), hold on, title('T')
subplot(3,1,2), plot(timead, rad2deg(alphaHS)), hold on, plot(timead, rad2deg(ualpha)), hold on, title('alpha')
subplot(3,1,3), plot(timead, rad2deg(betaHS)), hold on, plot(timead, rad2deg(ubeta)), hold on, title('beta')
figure()
plot(timead, rHS), hold on, plot(timead, r), hold on,
%%
function J = DTmethod(X, data)

    N =data.n_int;
    m = zeros(1, N);
    x = zeros(N,7);
    T = m;
        time = data.time;
% alpha = m; beta = m;


    for k=1:N
        x(k,:) = X((k-1)*10 + 1:(k-1)*10+7);
        m(k) = x(k,7);
        T(k) = X((k-1)*10+8); 
%         alpha(k) = X((k-1)*10+9);
%         beta(k) = X((k-1)*10+10);
    end
    
    J = 0;
    for k = 2:N-1
        %integration of cost function
        dt = time(k+1) - time(k);
        J = J +  dt/6 * (T(k-1) + 4*T(k) + T(k + 1));
    end
    J = J + T(end)*dt;
%     T_inplane = a_in .*m * 1000;
%     T_outplane = a_out .*m * 1000;
%     T = (T_inplane.^2 + T_outplane.^2).^0.5;
%     J = m(1) - m(end);

end

function [c, ceq] = EoM(X, data)
    muS = data.muS;
    time = data.time;
    N = data.n_int;
        
    x = zeros(N,7);
    m = zeros(1, N);
    T = m; ualpha = m; ubeta = m;
    
    for k=1:N
        %adimensional variables
        x(k,:) = X((k-1)*10 + 1:(k-1)*10+7);
        % dimensional variables
        T(k) = X((k-1)*10+8); 
        ualpha(k) = X((k-1)*10+9);
        ubeta(k) = X((k-1)*10+10); 
        
        m(k) = x(k,7);
    end
 
    %initial and final states
    xi = data.xi;
    xf = data.xf;
    BCi = xi(1:6) - x(1,1:6);
    BCf = xf(1:6) - x(end,1:6);
    
    Xprop = zeros(N, 7);
    Xprop(1,:) = x(1,:); %first value equal
    
    %RK4 backward integration (in each interval of delta t)
    for k = 1:N-1
%         xx = Xprop(k,:);
        xx = x(k,:);

        hhh = time(k+1) - time(k);
        K1 = polar_ad(time(k), xx, T(k), ualpha(k), ubeta(k), muS, data);
        K2 = polar_ad(time(k)+ hhh/2, xx + hhh/2*K1, 0.5*(T(k+1) + T(k)), 0.5*(ualpha(k+1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)), muS, data);
        K3 = polar_ad(time(k)+ hhh/2, xx + hhh/2*K2, 0.5*(T(k+1) + T(k)), 0.5*(ualpha(k+1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)),  muS, data);
        K4 = polar_ad(time(k)+ hhh, xx + hhh*K3, T(k+1), ualpha(k+1), ubeta(k+1), muS, data);    

        Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
%         Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end
%     for k = N:-1:2
%         xx = Xprop(k,:);
% %         xx = x(k,:);
% 
%         hhh = time(k) - time(k-1);
%         K1 = polar_ad(time(k), xx, T(k), ualpha(k), ubeta(k), muS, data);
%         K2 = polar_ad(time(k)+ hhh/2, xx + hhh/2*K1, 0.5*(T(k-1) + T(k)), 0.5*(ualpha(k-1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)), muS, data);
%         K3 = polar_ad(time(k)+ hhh/2, xx + hhh/2*K2, 0.5*(T(k-1) + T(k)), 0.5*(ualpha(k-1) + ualpha(k)), 0.5*(ubeta(k)+ubeta(k+1)),  muS, data);
%         K4 = polar_ad(time(k)+ hhh, xx + hhh*K3, T(k+1), ualpha(k-1), ubeta(k-1), muS, data);    
% 
% %         Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
%         Xprop(k-1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
%     end
    csi = zeros(1, (N-1)*7);
    %definition of the initial guess (sub-optimal conway solution)
    for ii = 1:N-1
        csi((ii-1)*7 +1 : (ii-1)*7 +7) = (x(ii+1,:) - Xprop(ii+1,:));
    end
  
    %collocation constraints
    DDD = zeros(1,(N-1)*7);
    
    %computation of collocation point
    for k = 1:N-1

        h = time(k+1) - time(k);

        xk = x(k,:); 
        xkk = x(k+1,:); 
        yk  = polar_ad(time(k), xk, T(k), ualpha(k), ubeta(k), muS, data);
        ykk  = polar_ad(time(k+1), xkk, T(k+1), ualpha(k+1), ubeta(k+1), muS, data);

        xc = 0.5*(xk + xkk) + h/8*(yk-ykk);    
        uc = ([T(k) ualpha(k) ubeta(k)] + [T(k+1) ualpha(k+1) ubeta(k+1)])/2;
        
        yc = polar_ad(time(k+1), xc, uc(1), uc(2), uc(3), muS, data);

        DDD((k-1)*7 + 1:(k-1)*7 + 7) =  xk - xkk + h/6*(yk + 4*yc + ykk);
    
    end   
    
    %quality constrataints
    ceq = [csi, DDD]; 
    c = [BCi BCf];
%     c = [];
end

function [lb, ub] = LBUB(XX0, data)
%X0 vector of propagated ODE45
    N = data.n_int;
    T_lb = 0;
    T_ub = data.Tmax;
    alpha_lb = 0;
    alpha_ub = 2*pi;
    beta_lb = -pi/2;
    beta_ub = pi/2;
    
    %limits definition
    ub = zeros(1,10*N);
    lb = zeros(1,10*N);
    
    %r
    r_ub = 2;
    r_lb = 0.8;

    %theta
    th_ub = Inf;
    th_lb = 0;
    %z
    z_ub =  0.05;
    z_lb =  -0.05;
    %vr
    vr_ub = 2;
    vr_lb = -2;
    %theta_dot
    thd_ub = Inf;
    thd_lb = -Inf;
    %vz
    vz_ub = 1;
    vz_lb = -1;
    
    %
    m_ub = XX0(7)*1.1;
    m_lb = 0.98;
    
    for k = 2:N-1
        ub((k-1)*10 +1) = r_ub;
        ub((k-1)*10 +2) = th_ub;
        ub((k-1)*10 +3) = z_ub;
        ub((k-1)*10 +4) = vr_ub;
        ub((k-1)*10 +5) = thd_ub;
        ub((k-1)*10 +6) = vz_ub;
        ub((k-1)*10 +7) = m_ub;
        ub((k-1)*10 +8) = T_ub;
        ub((k-1)*10 +9) = alpha_ub;
        ub((k-1)*10 +10) = beta_ub;
    
        lb((k-1)*10 +1) = r_lb;
        lb((k-1)*10 +2) = th_lb;
        lb((k-1)*10 +3) = z_lb;
        lb((k-1)*10 +4) = vr_lb;
        lb((k-1)*10 +5) = thd_lb;
        lb((k-1)*10 +6) = vz_lb;
        lb((k-1)*10 +7) = m_lb;%data.Mdry;
        lb((k-1)*10 +8) = T_lb;
        lb((k-1)*10 +9) = alpha_lb;
        lb((k-1)*10 +10) = beta_lb;
    end
    
    %Initial condition
    ub(1:10) = XX0(1:10); 
    ub(8) = T_ub; ub(9) = alpha_ub; ub(10) = beta_ub;
    ub(7) = m_ub;
    lb(1:10) = XX0(1:10); 
    lb(7) = m_lb;
    lb(8) = 0; lb(9) = alpha_lb; lb(10) = beta_lb;
    
    %Final condition
    ub(end-9:end) = XX0(end-9:end);
    ub(end-2)=T_ub; ub(end-1)=alpha_ub; ub(end)=beta_ub;
    ub(end-3) = m_ub;
    lb(end-9:end) = XX0(end-9:end);
    lb(end-2)=0; lb(end-1)=alpha_lb; lb(end)=beta_lb;
    lb(end-3) = m_lb;
end