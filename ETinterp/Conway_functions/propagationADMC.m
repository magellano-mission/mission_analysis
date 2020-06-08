function [Xd] = propagationADMC(T, ualpha, ubeta, X, muS, time, data)

N = data.n_int;
DU = astroConstants(2);
TU = (DU^3/muS).^0.5;
MU = data.Mdry;

Xad = zeros(N,7); Xd = Xad;
%adimensionalization of initial conditions
Xad(:,1) = X(:,1)/DU;
Xad(:,2) = X(:,2); %already adimensional
Xad(:,3) = X(:,3)/DU;
Xad(:,4) = X(:,4)/DU*TU; 
Xad(:,5) = X(:,5)*TU; 
Xad(:,6) = X(:,6)/DU*TU;
Xad(:,7) = X(:,7)/MU;

data.muS = muS;   
timead = time/TU;

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


end

