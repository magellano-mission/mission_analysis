function y = ga_conway(x, data_stacks)
t0 = x(1); 
TOF = x(2);
N_rev = round(x(3));
q = x(4);

data_stacks.t0sym = t0; %mjd2000 

[kepEarth, muS] = uplanet(t0, 3);
[rE0, vE0] = kep2car2(kepEarth, muS);
[kepMars, muS] = uplanet(t0 + TOF, 4);
[rMend, vMend] = kep2car2(kepMars, muS);

[~, ~, ~, ~, ~, ~, ~, ~, RCRRv, ~, vr, vt, v1perp, v2perp, ~, ~, ~, ~, ~, m, T, ~] = ...
    Conway(t0, TOF, N_rev, q, data_stacks);

if ~isnan(m)
r1vers = rE0/norm(rE0);
r2vers = rMend/norm(rMend) ;
dirt1 = cross(RCRRv , r1vers);
dirt2 = cross(RCRRv , r2vers);

vr1 = r1vers * vr(1);
vr2 = r2vers * vr(end);

vt1 = dirt1 * vt(1);
vt2 = dirt2 * vt(end);

v1 = vr1 + vt1 + v1perp; v2 = vr2 + vt2 + v2perp;
v_inf1 = v1 - vE0;
v_inf2 = vMend - v2;

y(1) = max(abs(T));
y(2) = m(1) - m(end);
y(3) = norm(v_inf2);
y(4) = norm(v_inf1);
else
    y = 999999*(1*rand(4,1));
end
end

