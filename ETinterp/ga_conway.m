function y = ga_conway(x, data_stacks)
t0 = x(1);              TOF = x(2);
N_rev = round(x(3));      q = x(4);
v_inf = x(5); alpha = x(6); beta = x(7);

[kepEarth, muS] = uplanet(t0      ,3);
[kepMars, ~]    = uplanet(t0 + TOF,4);

[RE, vE] = kep2car2(kepEarth, muS); %km....
[R2, v2] = kep2car2(kepMars, muS);  %km....

R1 = RE;

%definition of plane of motion
r1norm = norm(R1);
r2norm = norm(R2);

r1vers = R1/r1norm;
r2vers = R2/r2norm;

RIvcRFv = cross(r1vers, r2vers);
RCRRv = RIvcRFv/norm(RIvcRFv);

if RIvcRFv(3) <0
    RCRRv = -RCRRv;
end

%adding TMI maneuver
v1 = vE + v_inf*(sin(beta)*cos(alpha)*r1vers + ...
                 sin(beta)*sin(alpha)*cross(RCRRv,r1vers) + ...
                 cos(beta)*RCRRv);

[ m, T ] = Conway(TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, RCRRv, RIvcRFv, v1, v2, muS, data_stacks);

if ~isnan(m)

y(1) = max(abs(T));
y(2) = m(1) - m(end);

else
    y = 999999*(1*rand(2,1));
end
end

