function t_l = MarsETphasideal(tend, r0, m, T, mi, dth_phas)

a = T/m;
n0 = sqrt(mi/r0^3);
t = 1:tend;

rt = r0./(1 - a*t*sqrt(r0/mi));
rf = rt(end);
nt = sqrt(mi./rt.^3);
intf = nt - n0;
dth = trapz(t, intf);

th_l = dth_phas - 2*dth;
nf = sqrt(mi/rf^3);
t_l = th_l/(nf - n0);