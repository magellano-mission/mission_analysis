function [X] = EoMpropRK4(X0, time, a_in, a_out, muS, data, plotsflag, r, TH, z, vr, th_dot, vz, m)
%propagation of function
    if nargin == 6
       plotsflag = 0; 
    end
    X = zeros(data.n_int, 7);
    X(1,:) = X0;
    %RK4 forward integration
    for k = 1:data.n_int-1
        x = X(k,:);
        
        hhh = time(k+1) - time(k);
        K1 = EoMpolar(time(k), x, a_in(k), a_out(k), muS, data);
        K2 = EoMpolar(time(k)+ hhh/2, x + hhh/2*K1, 0.5*(a_in(k+1) + a_in(k)), 0.5*(a_out(k+1) + a_out(k)), muS, data);
        K3 = EoMpolar(time(k)+ hhh/2, x + hhh/2*K2, 0.5*(a_in(k+1) + a_in(k)), 0.5*(a_out(k+1) + a_out(k)), muS, data);
        K4 = EoMpolar(time(k)+ hhh, x + hhh*K3, a_in(k+1), a_out(k+1), muS, data);    
        
        X(k+1,:) = X(k,:) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end

    if plotsflag
        TX_inplane = a_in .*X(:,7)' * 1000;
        TX_outplane = a_out .*X(:,7)' * 1000;
        TX = (TX_inplane.^2 + TX_outplane.^2).^0.5;
        figure()
        sgtitle('States validation - Dynamics states')
        subplot(6,2,1), plot(TOFr, r,'k:'), ylabel('r [km]'), title('Conway Solution')
        subplot(6,2,2), plot(TOFr, X(:,1),'k:'), ylabel('r [km]'), title('Propagated Solution')
        subplot(6,2,3), plot(TOFr, TH, 'k:'), ylabel('$\theta$ [km]'), 
        subplot(6,2,4), plot(TOFr, X(:,2),'k:'), ylabel('$\theta$ [km]')
        subplot(6,2,5), plot(TOFr, z, 'k:'),  ylabel('z[km]')
        subplot(6,2,6), plot(TOFr, X(:,3),'k:'), ylabel('z[km]')
        subplot(6,2,7), plot(TOFr, vr, 'k:'), ylabel('$v_r$ [km]')
        subplot(6,2,8), plot(TOFr, X(:,4), 'k:'), ylabel('$v_r$ [km]')
        subplot(6,2,9), plot(TOFr, r.*th_dot, 'k:'),  ylabel('$\dot{\theta}$ [km]')
        subplot(6,2,10), plot(TOFr, X(:,1).*X(:,5), 'k:'), ylabel('$\dot{\theta}$[km]')
        subplot(6,2,11), plot(TOFr, vz, 'k:'),  xlabel('TOFr [days]'), ylabel('$v_z$ [km]') 
        subplot(6,2,12), plot(TOFr, X(:,6), 'k:'),xlabel('TOFr [days]'), ylabel('$v_z$ [km]')

        figure(), sgtitle('States validation - Mass and Thrust')
        subplot(2,2,1), plot(TOFr, m, 'k:'), xlabel('TOF [days]'), ylabel('Spacecraft Mass [kg]')
        title('Conway Solution'), 
        subplot(2,2,2), plot(TOFr, X(:,7), 'k:'), xlabel('TOF [days]'), ylabel('Spacecraft Mass [kg]')
        title('Propagated Solution')
        subplot(2,2,3), plot(TOFr, T, 'k:'),  xlabel('TOF [days]'), ylabel('Thrust [N]')
        subplot(2,2,4), plot(TOFr, TX, 'k:'),  xlabel('TOF [days]'), ylabel('Thrust [N]')
    end
end

