function [dt, Jmin] = ETphasOptim(T, lb, ub, data)
%{
Function to get the best Electric phasing maneuver
The output dt is the total period of the maneuver in days.
The total period is perfectly divided in 2 segment of continuos thrust;
no free motion between those segment is added, in order to get the minimum tof.
%}

NN = 5;
j = 0;
Jmin = 1e5;
Jold = zeros(NN, 1);
tol = 1e-3;
r = NN + 1;

while Jmin > tol && r >= NN
    J = Jold;
    j = j+1;
    if j~=1
        for i = 2:NN-1
            TT = linspace(lb, ub, NN);
            data.T = T;
            [~, Y1] = ode113(@ETphaseIntegration, [0, TT(i)], data.Y0, data.opt, data);
            % chaser dynamics
            data.T = -T;
            [~, Y3] = ode113(@ETphaseIntegration, [0, TT(i)], Y1(end, :), data.opt, data);
            
            % target dynamics
            [~, Y4] = ode113(@ETphaseIntegration, [0, 2*TT(i)], data.Yfree, data.opt, data);
            
            J(i) = norm(Y3(end, 1:3) - Y4(end, 1:3));
        end
        
    else
        for i = 1:NN
            TT = linspace(lb, ub, NN);
            data.T = T;
            % chaser dynamics
            [~, Y1] = ode113(@ETphaseIntegration, [0, TT(i)], data.Y0, data.opt, data);
            data.T = -T;
            [~, Y3] = ode113(@ETphaseIntegration, [0, TT(i)], Y1(end, :), data.opt, data);
            
            % target dynamics
            [~, Y4] = ode113(@ETphaseIntegration, [0, 2*TT(i)], data.Yfree, data.opt, data);
            
            J(i) = norm(Y3(end, 1:3) - Y4(end, 1:3));
        end
    end
    
    [Jmin, ind] = min(J);
    
    if ind ~= 1
        lb = floor(TT(ind-1));
        Jlb = J(ind-1);
    else
        Jlb = J(1);
    end
    
    if ind ~= NN
        ub = ceil(TT(ind+1));
        Jub = J(ind+1);
    else
        Jub = J(NN);
    end
    
    Jold = [Jlb; 0; 0; 0; Jub];
    
    r = (ub - lb);
end

dt = 2*TT(ind)/86400;
