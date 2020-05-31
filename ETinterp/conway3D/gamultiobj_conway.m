function y = gamultiobj_conway(x, data_stacks)
%variables
% t0 = x(1);              TOF = x(2);
% N_rev = round(x(3));      q = x(4);
% v_inf = x(5); alpha = x(6); beta = x(7);

%variables
% %arrival at the date of the first plus 6 months
% tf = 9.619974167516721e+03 + 8.361779091278747e+02 + 360;  TOF = x(1);
% t0 = tf - TOF;
% N_rev = round(x(2));      q = x(3);
% v_inf = x(4); alpha = x(5); beta = x(6);

%variables
% N_rev = round(x(1));      q = x(2);
% v_inf      =          1.019932213771344;
% alpha      =          1.339209562214939;
% beta       =          2.689364771171173;
% t0 = 9.618881860211146e+03 + 6; TOF = 9.084518115422542e+02;

[kepEarth, muS] = uplanet(t0,3);
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

%vinf @ Mars = 0;
             
[ m, T ] = Conway(TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, hvers, hh, v1, v2, muS, data_stacks);

if ~isnan(m) 

        y(1) = max(abs(T))*10000;
        y(2) = m(1) - m(end);
else
    y = 1e6*ones(2,1);
end
end

