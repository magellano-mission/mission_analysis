function [J, OmF, iFout, DOmout] = BCforCapturing(x, orbit)

Om0 = orbit.Om0;
iF = orbit.iF;
DOm = orbit.DOm;
om0 = orbit.om;
% opt = orbit.opt;

th0 = x(1);
i0 = x(2);
Dt = x(3);

[~, Y] = ode45(@IntegrGauss, [0 Dt], [i0, Om0, th0, om0], [], orbit);

iFout = Y(end, 1);
OmF = Y(end, 2);

DOmout = OmF - Om0;


J = abs(DOm - DOmout) + abs(iF - iFout);




