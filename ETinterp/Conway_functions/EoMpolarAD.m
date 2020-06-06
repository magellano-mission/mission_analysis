function [y] = EoMpolarAD(t, x, T, alpha, beta, muS, data)
r = x(1); th = x(2); z = x(3); vr = x(4); thd = x(5); vz = x(6); m = x(7);
s = (r^2 + z^2)^0.5;
DU = 149597870.7;
TU = (DU^3/muS).^0.5;
MU = data.Mdry;
muS = 1;

% DU = 1; 
% TU = 1;
% MU = 1;

m = m*MU;

y(1) = vr;
y(2) = thd;
y(3) = vz;
y(4) = -muS/s^3 * r + r*thd^2 + T/1000/m/DU*TU^2 *cos(beta)*sin(alpha);
y(5) = 1/r*(T/1000/m/DU*TU^2*cos(beta)*cos(alpha) - 2*vr*thd);
y(6) = -muS/(s)^(3) * z + T/1000/m/DU*TU^2*sin(beta);
y(7) = (-T/(data.Isp(T)*9.81))/MU*TU; 

end

