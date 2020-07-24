function [J, Y, T, VectorThrust, VectorLight, VectorThRange] = BCfindRHyp(x, data)

rp = x(1);                                % pericenter of arrival hyp [km]

e = 1 + rp*data.Vinf^2/data.mi;           % eccentricity of the hyp
p = rp*(e + 1);                           % semilatum rectus [km]
data.th0 = -acos((p/data.R_SOI - 1)/e);   % theta infinite

% Initial conditions
data.e_hyp = e;                           % eccentricity
data.sma = p/(data.e_hyp^2 - 1);          % semi-major axis [km] 

Kep0 = [data.sma, data.e_hyp, data.i, data.Om0, 0, data.th0];     
[R0, V0] = kep2car(Kep0, data.mi);        % cartesian IC
M0 = data.M;
Y0 = [R0; V0; M0];

[T, Y] = ode113(@ETcaptIntegration, [0, 1e8], Y0, data.opt, data);

VectorThrust = zeros(length(T),1);
VectorLight = zeros(length(T),1);
VectorThRange = zeros(length(T),1);
for k = 1:length(T)
    [~, parout] = ETcaptIntegration(T(k),Y(k,:)', data);
    VectorThrust(k) = parout{1};
    VectorLight(k) = parout{2};
    VectorThRange(k) = parout{3};
end

rFout = norm(Y(end,1:3));

J = rFout - data.r_final;

end