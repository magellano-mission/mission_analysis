%direct transcription

options = optimoptions('fmincon','Display','iter', ...
                        'MaxFunctionEvaluations',1e4,'Algorithm','interior-point');
XX0 = zeros(1, length(TH)*10);

%definition of the initial guess (sub-optimal conway solution)
for ii = 1:length(TH)
    XX0((ii-1)*10 +1: (ii-1)*10+10) = [X(ii,:) acc_inplane(ii) acc_out(ii) gamma(ii)];
end

[XSOL,fval,exitflag,output] = fmincon(@DCmethod, XX0,[],[],[],[],[],[],@mycon, options, data_stacks);

XHS = zeros(data_stacks.n_int, 7);
acc_inplaneHS = zeros(1,data_stacks.n_int);
acc_outHS=acc_inplaneHS;
gammaHS=acc_inplaneHS;
%reconstruction of states
for ii = 1:length(TH)
    [XHS(ii,:)] = XSOL((ii-1)*10 +1: (ii-1)*10+7);
    [acc_inplaneHS(ii)] = XSOL((ii-1)*10 +8);
    [acc_outHS(ii)] = XSOL((ii-1)*10+9);
    [gammaHS(ii)] = XSOL((ii-1)*10+10);
end
%%
rHS = XHS(:,1)';
mHS = XHS(:,7)';
T_inplaneHS = acc_inplaneHS .*mHS * 1000;
T_outplaneHS = acc_outHS .*mHS * 1000;
THS = (T_inplaneHS.^2 + T_outplaneHS.^2).^0.5;
%%
figure()
sgtitle('Control allocation')
subplot(3,1,1), plot(time, acc_inplane), hold on, plot(time, acc_inplaneHS), hold on,
subplot(3,1,2), plot(time, acc_out), hold on, plot(time, acc_outHS), hold on,
subplot(3,1,3), plot(time, THS), hold on, plot(time, T), hold on,
figure()
plot(time, rHS), hold on, plot(time, r), hold on,
%%
function J = DCmethod(X, data)

    N =data.n_int;
    m = zeros(1, N);
    x = zeros(N,7);
    a_in = m; a_out = m; gamma = m;
    for k=1:N
    x(k,:) = X((k-1)*10 + 1:(k-1)*10+7);
    m(k) = x(k,7);
    a_in(k) = X((k-1)*10+8); a_out(k) = X((k-1)*10+9); gamma(k) = X((k-1)*10+10);
    end

%     T_inplane = a_in .*m * 1000;
%     T_outplane = a_out .*m * 1000;
%     T = (T_inplane.^2 + T_outplane.^2).^0.5;
%     J = max(abs(T));
    J = m(1) - m(end);

end

function [c, ceq] = mycon(X, data)
    muS = data.muS;
    time = data.time;
    N =data.n_int;
    m = zeros(1, N);
    x = zeros(N,7);
    a_in = m; a_out = m; gamma = m;
    for k=1:N
        x(k,:) = X((k-1)*10 + 1:(k-1)*10+7);
        m(k) = x(k,7);
        a_in(k) = X((k-1)*10+8); a_out(k) = X((k-1)*10+9); gamma(k) = X((k-1)*10+10);
    end
    
    %initial and final states
    xi = data.xi;
    xf = data.xf;
    BCi = xi - x(1,:);
    BCf = xf - x(end,:);
    
    T_inplane = a_in .*m * 1000;
    T_outplane = a_out .*m * 1000;
    T = (T_inplane.^2 + T_outplane.^2).^0.5;
      
    Xprop = zeros(N, 7);
    Xprop(1,:) = x(1,:); %first value equal
    %RK4 forward integration (in each interval of delta t)
    for k = 1:N-1
%         xx = Xprop(k,:);
        xx = x(k,:);

        hhh = time(k+1) - time(k);
        K1 = propagate_conway(time(k), xx, a_in(k), a_out(k), gamma(k), muS, data);
        K2 = propagate_conway(time(k)+ hhh/2, xx + hhh/2*K1, 0.5*(a_in(k+1) + a_in(k)), 0.5*(a_out(k+1) + a_out(k)), 0.5*(gamma(k+1)-gamma(k)), muS, data);
        K3 = propagate_conway(time(k)+ hhh/2, xx + hhh/2*K2, 0.5*(a_in(k+1) + a_in(k)), 0.5*(a_out(k+1) + a_out(k)), 0.5*(gamma(k+1)-gamma(k)), muS, data);
        K4 = propagate_conway(time(k)+ hhh, xx + hhh*K3, a_in(k+1), a_out(k+1), gamma(k+1), muS, data);    

%         Xprop(k+1,:) = Xprop(k,:) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
        Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end

    csi = zeros(1, (N-1)*7);
    %definition of the initial guess (sub-optimal conway solution)
    for ii = 1:N-1
        csi((ii-1)*7 +1 : (ii-1)*7+7) = (x(ii+1,:) - Xprop(ii+1,:));
    end

    
    %collocation constraints
    DDD = zeros(1,(N-1)*7);

    for k = 1:N-1
    
    xk = x(k,:); 
    a_ink  = a_in(k); a_outk = a_out(k); gammak = gamma(k);
    
    xkk = x(k+1,:); 
    a_inkk = a_in(k+1);     a_outkk = a_out(k+1);  gammakk = gamma(k+1);
    
    h = time(k+1) - time(k);
        
    yk  = propagate_conway(time(k), xk, a_ink, a_outk, gammak, muS, data);
    ykk = propagate_conway(time(k+1), xkk, a_inkk, a_outkk, gammakk, muS, data);
    
    %computation of collocation point
    xc = 0.5*(xk + xkk) + h/8*(yk-ykk);
    
    uc = ([a_ink a_outk gammak] + [a_inkk a_outkk gammakk])/2;

    yc = propagate_conway(time(k), xc, uc(1), uc(2), uc(3), muS, data);

    DDD((k-1)*7 + 1:(k-1)*7 + 7) =  xk - xkk + h/6*(yk + 4*yc + ykk);
    
    end   
    
    c = [];
    ceqT = zeros(1,3*N);
    for cc = 1:N
        if abs(T(cc)) >= data.Tmax %thrust on
            ceqT((cc-1)*cc+1:(cc-1)*cc+3) = [(abs(T(cc)) - data.Tmax), ...
                                            (abs(a_in(cc))-(data.Tmax/1000/m(cc))), ...
                                            (abs(a_out(cc))-(0/1000/m(cc)))];
        else
            ceqT((cc-1)*cc+1:(cc-1)*cc+3) = [abs(T(cc)), ...
                                             abs(a_in(cc)), ...
                                             abs(a_out(cc))];
        end
    end
    
    
    %quality constrataints
    ceq = [csi, BCi, BCf, DDD, ceqT]; 

end