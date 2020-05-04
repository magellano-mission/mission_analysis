function Y = trajectory(X0, data)

Y = zeros(data.NT, 3);
tspan = data.tspan;
th0 = X0(6);
a0 = X0(1);
e0 = X0(2);
i0 = X0(3);


for kk = 1:data.NT
    X0(1) = a0 + data.adot*(tspan(kk) - tspan(1));
    X0(2) = e0 + data.edot*(tspan(kk) - tspan(1));
    X0(3) = i0 + data.idot*(tspan(kk) - tspan(1));
    
    M = th0 + data.n*(tspan(kk) - tspan(1));
    if data.edot == 0
        X0(6) = M;
    else
        X0(6) = M + (2*X0(2) - 1/4*X0(2)^2)*sin(M) + 5/4*X0(2)^2*sin(2*M) + 13/12*X0(2)^3*sin(3*M);
    end
        
    Y(kk, :) = kep2car(X0, data.mi);
end


end