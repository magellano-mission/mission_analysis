function [T, Y] = trajectory(X0, tspan, mi, opt, a, perturb)

NT = length(tspan);
Y = zeros(NT, 3);
if perturb
    [~, X] = ode113(@gauss, tspan, X0, opt);

    for kk = 1 : NT
        [rr, ~] = kep2car(X(kk,:), mi);
        Y(kk,:) = rr';
    end
else
    th0 = X0(6);
    n = sqrt(mi/a^3);
    for kk = 1:NT
        X0(6) = th0 + n*(tspan(kk) - tspan(1));
        Y(kk, :) = kep2car(X0, mi);
    end
    
end

T = tspan;


end