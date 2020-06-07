function [J, rFout] = BCfindRHyp(x, data)

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

[~, Y] = ode113(@ETcaptIntegration, [0, data.TT], Y0, data.opt, data);

rFout = norm(Y(end,1:3));

J = abs(rFout - data.r_final);

end