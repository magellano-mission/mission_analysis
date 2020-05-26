%direct transcription

options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
XX0 = zeros(1, length(TH)*10);
%definition of the initial guess (sub-optimal conway solution)
for ii = 1:length(TH)
    XX0((ii-1)*10 +1: (ii-1)*10+10) = [X(ii,:) acc_inplane(ii) acc_out(ii) gamma(ii)];
end
x = fmincon(@DCmethod, XX0,[],[],[],[],[],[],@mycon, options, data_stacks);

% XHS = zeros(data_stacks.n_int, 7);
% acc_inplaneHS = zeros(1,data_stacks.n_int);
% acc_outHS=acc_inplaneHS;
% gammaHS=acc_inplaneHS;
% %reconstruction of states
% for ii = 1:length(TH)
%     [XHS(ii,:), acc_inplaneHS(ii), acc_outHS(ii), gammaHS(ii)] = XSOL((ii-1)*10 +1: (ii-1)*10+10);
% end
%%

function D = DCmethod(X, data)
    muS = data.muS;
    N = data.n_int;
    time = data.time;
    DDD = zeros(1,N*7);
    
    for k = 1:N-1
    
    xk = X((k-1)*10 + 1:(k-1)*10+7); 
    acc_inplanek  = X((k-1)*10+8); acc_outk = X((k-1)*10+9); gammak = X((k-1)*10+10);
    
    xkk = X(k*10+1:k*10+7); 
    acc_inplanekk = X(k*10+8);     acc_outkk = X(k*10+9);    gammakk = X(k*10+10);
    
    h = time(k+1) - time(k); %x fixed, u = [acc_in acc_out gamma]
        
    yk  = propagate_conway(time(k), xk, acc_inplanek, acc_outk, gammak, muS, data);
    ykk = propagate_conway(time(k+1), xkk, acc_inplanekk, acc_outkk, gammakk, muS, data);
    
    %computation of collocation point
    xc = 0.5*(xk + xkk) + h/8*(yk-ykk);
    
    uc = ([acc_inplanek acc_outk gammak] + [acc_inplanekk acc_outkk gammakk])/2;

    yc = propagate_conway(time(k), xc, uc(1), uc(2), uc(3), muS, data);
    
    
    DDD((k-1)*7 + 1:(k-1)*7 + 7) = [ xk - xkk + h/6*(yk + 4*yc + ykk)];
    
    end    
    D = norm(DDD);
end


function [c, ceq] = mycon( X, data)
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

    T_inplane = a_in .*m * 1000;
    T_outplane = a_out .*m * 1000;
    T = (T_inplane.^2 + T_outplane.^2).^0.5;
    c = abs(T) - data.Tmax;

    Xprop = zeros(N, 7);
    Xprop(1,:) = x(1,:); %first value equal
    %RK4 forward integration
    for k = 1:N-1
        xx = Xprop(k,:);

        hhh = time(k+1) - time(k);
        K1 = propagate_conway(time(k), xx, a_in(k), a_out(k), gamma(k), muS, data);
        K2 = propagate_conway(time(k)+ hhh/2, xx + hhh/2*K1, 0.5*(a_in(k+1) + a_in(k)), 0.5*(a_out(k+1) + a_out(k)), 0.5*(gamma(k+1)-gamma(k)), muS, data);
        K3 = propagate_conway(time(k)+ hhh/2, xx + hhh/2*K2, 0.5*(a_in(k+1) + a_in(k)), 0.5*(a_out(k+1) + a_out(k)), 0.5*(gamma(k+1)-gamma(k)), muS, data);
        K4 = propagate_conway(time(k)+ hhh, xx + hhh*K3, a_in(k+1), a_out(k+1), gamma(k+1), muS, data);    

        Xprop(k+1,:) = Xprop(k,:) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end

    csi = zeros(1, N*7);
    %definition of the initial guess (sub-optimal conway solution)
    for ii = 1:N
        csi((ii-1)*7 +1 : (ii-1)*7+7) = (x(ii,:)-Xprop(ii,:));
    end
    %initial and final states
    xi = data.xi;
    xf = data.xi;
    BCi = xi - x(1,:);
    BCf = xf - x(end,:);

    ceq = [csi, BCi, BCf ]; %add propagation of EoM

end