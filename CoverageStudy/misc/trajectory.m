function [T, Y] = trajectory(X0, tspan, data, opt, curr)

Y = zeros(data.NT, 3);
if data.perturb
    [~, X] = ode113(@gauss, tspan, X0, opt);

    for kk = 1 : data.NT
        [rr, ~] = kep2car(X(kk,:), data.mi);
        Y(kk,:) = rr';
    end
else
    th0 = X0(6);
    n = sqrt(data.mi/curr.sma^3);             % angular velocity
    for kk = 1:data.NT
        X0(6) = th0 + n*(tspan(kk) - tspan(1));
        Y(kk, :) = kep2car(X0, data.mi);
    end
    
end

T = tspan;


end