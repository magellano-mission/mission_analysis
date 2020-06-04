function [y] = EoMpolar2(t, x, T, alpha, beta, muS, data)
r = x(1); th = x(2); z = x(3); vr = x(4); thd = x(5); vz = x(6); m = x(7);
s = (r^2 + z^2)^0.5;

y(1) = vr;
y(2) = thd;
y(3) = vz;
y(4) = -muS/s^3 * r + r*thd^2 + T/1000/m*cos(beta)*sin(alpha);
y(5) = 1/r*(T/1000/m*cos(beta)*cos(alpha) - 2*vr*thd);
y(6) = -muS/(s)^(3) * z + T/1000/m*sin(beta);
y(7) = (-T/(data.Isp*9.81)); 

end

