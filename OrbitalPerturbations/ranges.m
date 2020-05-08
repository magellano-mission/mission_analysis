function [kep_ranges] = ranges(state, h_range, n, n_orb)

mu_m = astroConstants(14);
R_m = 3389.5;

hh = linspace(h_range(1), h_range(2), n);

kep_ranges = zeros(n,4);

for kk = 1 : n
    

a = R_m + hh(kk);
state(1) = a;

T = 2 * pi * sqrt(a^3 / mu_m);
t_span = [0 n_orb*T];

[~, x] = ode113(@gauss, t_span,  state);


kep_min = min(x, [], 1);
kep_max = max(x, [], 1);

kep_ranges(kk,:) = [kep_min(3), kep_max(3), kep_min(4), kep_max(4)];
end

end