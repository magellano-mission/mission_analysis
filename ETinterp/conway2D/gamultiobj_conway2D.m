function y = gamultiobj_conway2D(x, data_stacks)
t0 = x(1);
TOF = x(2);
N_rev = round(x(3));  
v_inf = x(4); alpha = x(5);
% v_infcap = x(8); alphacap = x(9); betacap = x(10);
beta = pi/2;
[kepEarth, muS] = uplanet(t0      ,3);
[kepMars, ~]    = uplanet(t0 + TOF,4);

[R1, v1] = kep2car2(kepEarth, muS); %km....
[R2, v2] = kep2car2(kepMars, muS);  %km....

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
v1 = v1 + v_inf*(sin(beta)*cos(alpha)*r1vers + ...
                 sin(beta)*sin(alpha)*cross(RCRRv,r1vers) + ...
                 cos(beta)*RCRRv);

% %adding v_inf at Mars capture
% v2 = v2 + v_infcap*(sin(betacap)*cos(alphacap)*r2vers + ...
%                  sin(betacap)*sin(alphacap)*cross(RCRRv,r2vers) + ...
%                  cos(betacap)*RCRRv);
             
[ m, T ] = Conway2D(TOF, N_rev, r1norm, r2norm, r1vers, r2vers, RCRRv, RIvcRFv, v1, v2, muS, data_stacks);

if ~isnan(m) 
 if max(abs(T)) < 0.25
        y(1) = 1e4*max(abs(T));
        y(2) = m(1) - m(end);
 else
     y = 1e6*ones(2,1);
 end
else
    y = 1e6*ones(2,1);
end
end

