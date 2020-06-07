function [X] = EoMpropRK4_2(X0, time, T, alpha, beta, muS, data, plotsflag, r, TH, z, vr, th_dot, vz, m, TOFr)
%propagation of function
    if nargin == 7
       plotsflag = 0; 
    end
    X = zeros(data.n_int, 7);
    X(1,:) = X0;
    N = data.n_int;
    %RK4 forward integration
    for k = 1:N-1
        xx = X(k,:);
%         xx = x(k,:);

        hhh = time(k+1) - time(k);
        K1 = EoMpolar2(time(k), xx, T(k), alpha(k), beta(k), muS, data);
        K2 = EoMpolar2(time(k)+ hhh/2, xx + hhh/2*K1, 0.5*(T(k+1) + T(k)), 0.5*(alpha(k+1) + alpha(k)), 0.5*(beta(k)+beta(k+1)), muS, data);
        K3 = EoMpolar2(time(k)+ hhh/2, xx + hhh/2*K2, 0.5*(T(k+1) + T(k)), 0.5*(alpha(k+1) + alpha(k)), 0.5*(beta(k)+beta(k+1)),  muS, data);
        K4 = EoMpolar2(time(k)+ hhh, xx + hhh*K3, T(k+1), alpha(k+1), beta(k+1), muS, data);    

%         Xprop(k+1,:) = xx + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
        X(k+1,:) = X(k,:) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end

    if plotsflag
        TX_inplane = a_in .*X(:,7)' * 1000;
        TX_outplane = a_out .*X(:,7)' * 1000;
        TX = (TX_inplane.^2 + TX_outplane.^2).^0.5;
        figure()
        sgtitle('States validation - Dynamics states')
        subplot(6,1,1), plot(TOFr, r'- X(:,1),'k:'), ylabel('r [km]'), title('Propagated Solution')
        subplot(6,1,2), plot(TOFr, TH'- X(:,2),'k:'), ylabel('$\theta$ [km]')
        subplot(6,1,3), plot(TOFr, z'- X(:,3),'k:'), ylabel('z[km]')
        subplot(6,1,4), plot(TOFr, vr'- X(:,4), 'k:'), ylabel('$v_r$ [km]')
        subplot(6,1,5), plot(TOFr, (r.*th_dot)'- X(:,1).*X(:,5), 'k:'), ylabel('$\dot{\theta}$[km]')
        subplot(6,1,6), plot(TOFr, vz'- X(:,6), 'k:'),xlabel('TOFr [days]'), ylabel('$v_z$ [km]')

        figure(), sgtitle('States validation - Mass and Thrust')
        subplot(2,2,1), plot(TOFr, m, 'k:'), xlabel('TOF [days]'), ylabel('Spacecraft Mass [kg]')
        title('Conway Solution'), 
        subplot(2,2,2), plot(TOFr, X(:,7), 'k:'), xlabel('TOF [days]'), ylabel('Spacecraft Mass [kg]')
        title('Propagated Solution')
        subplot(2,2,3), plot(TOFr, T, 'k:'),  xlabel('TOF [days]'), ylabel('Thrust [N]')
        subplot(2,2,4), plot(TOFr, TX, 'k:'),  xlabel('TOF [days]'), ylabel('Thrust [N]')
    end
end

