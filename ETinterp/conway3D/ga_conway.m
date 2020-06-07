function y = ga_conway(x, data_stacks)
t0 = x(1);              TOF = x(2);
N_rev = round(x(3));      q = x(4);
v_inf = x(5); alpha = x(6); beta = x(7);
v_infcap = x(8); alphacap = x(9); betacap = x(10);

[kepEarth, muS] = uplanet(t0      ,3);
[kepMars, ~]    = uplanet(t0 + TOF,4);

[R1, v1] = kep2car2(kepEarth, muS); %km....
[R2, v2] = kep2car2(kepMars, muS);  %km....

%definition of plane of motion
r1norm = norm(R1);
r2norm = norm(R2);

r1vers = R1/r1norm;
r2vers = R2/r2norm;

hh = cross(r1vers, r2vers);
hvers = hh/norm(hh);

if hh(3) <0
    hvers = -hvers;
end

%adding TMI maneuver
v1 = v1 + v_inf*(sin(beta)*cos(alpha)*r1vers + ...
                 sin(beta)*sin(alpha)*cross(hvers,r1vers) + ...
                 cos(beta)*hvers);

%adding v_inf at Mars capture
v2 = v2 + v_infcap*(sin(betacap)*cos(alphacap)*r2vers + ...
                 sin(betacap)*sin(alphacap)*cross(hvers,r2vers) + ...
                 cos(betacap)*hvers);
             
[ m, T ] = Conway(TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, hvers, hh, v1, v2, muS, data_stacks);

if ~isnan(m) 

        y(1) = max(abs(T));
        
else
    y = 1e5*ones(1,1);
end
end

