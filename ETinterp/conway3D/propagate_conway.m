function [y] = propagate_conway(t, x, a_in, a_out, gamma, muS, data)
r = x(1); th = x(2); z = x(3); vr = x(4); thd = x(5); vz = x(6); m = x(7);

T_inplane = a_in .*m * 1000;
T_outplane = a_out .*m * 1000;
T = (T_inplane.^2 + T_outplane.^2).^0.5;

y(1) = vr;
y(2) = thd;
y(3) = vz;
y(4) = -muS/r^2 + r*thd^2 + a_in*sin(gamma);
y(5) = 1/r*(a_in*cos(gamma) - 2*vr*thd);
y(6) = -muS/(r)^(3) * z + a_out;
y(7) = -abs(T)/(data.Isp*9.81);
end

