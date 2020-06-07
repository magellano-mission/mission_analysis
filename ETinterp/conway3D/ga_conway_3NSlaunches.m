function [y, T1, T2, T3, m1, m2, m3] = ga_conway_3NSlaunches(x, data_stacks)
%%% variables
%first launch
t01 = x(1);              TOF1 = x(2);
N_rev1 = (x(3));      q1 = x(4);
v_inf1 = x(5); alpha1 = x(6); beta1 = x(7);

%second launch
TOF2 = x(8);
N_rev2 = (x(9));      q2 = x(10);
v_inf2 = x(11); alpha2 = x(12); beta2 = x(13);
t02 = t01 + TOF1 + 186 - TOF2 + x(20);

%third launch
TOF3 = x(14);
N_rev3 = (x(15));      q3 = x(16);
v_inf3 = x(17); alpha3 = x(18); beta3 = x(19);
t03 = t01 + TOF1 + 372 - TOF3 + x(21);

[kepEarth_1, muS] = uplanet(t01,3);
[kepMars_1, ~]    = uplanet(t01 + TOF1,4);

[R1_1, v1_1] = kep2car2(kepEarth_1, muS); %km....
[R2_1, v2_1] = kep2car2(kepMars_1, muS);  %km....

%definition of plane of motion
r1norm_1 = norm(R1_1);
r2norm_1 = norm(R2_1);

r1vers_1 = R1_1/r1norm_1;
r2vers_1 = R2_1/r2norm_1;

hh_1 = cross(r1vers_1, r2vers_1);
hvers_1 = hh_1/norm(hh_1);

if hh_1(3) <0
    hvers_1 = -hvers_1;
end

%adding TMI maneuver
v1_1 = v1_1 + v_inf1*(sin(beta1)*cos(alpha1)*r1vers_1 + ...
                 sin(beta1)*sin(alpha1)*cross(hvers_1,r1vers_1) + ...
                 cos(beta1)*hvers_1);

%vinf @ Mars = 0;
             
[ m1, T1 ] = Conway(TOF1, N_rev1, q1, r1norm_1, r2norm_1, r1vers_1, r2vers_1, hvers_1, hh_1, v1_1, v2_1, muS, data_stacks);


%second launch
[kepEarth_2, ~]   = uplanet(t02       ,3);
[kepMars_2, ~]    = uplanet(t02 + TOF2,4);

[R1_2, v1_2] = kep2car2(kepEarth_2, muS); %km....
[R2_2, v2_2] = kep2car2(kepMars_2, muS);  %km....

%definition of plane of motion
r1norm_2 = norm(R1_2);
r2norm_2 = norm(R2_2);

r1vers_2 = R1_2/r1norm_2;
r2vers_2 = R2_2/r2norm_2;

hh_2 = cross(r1vers_2, r2vers_2);
hvers_2 = hh_2/norm(hh_2);

if hh_2(3) <0
    hvers_2 = -hvers_2;
end

%adding TMI maneuver
v1_2 = v1_2 + v_inf2*(sin(beta2)*cos(alpha2)*r1vers_2 + ...
                 sin(beta2)*sin(alpha2)*cross(hvers_2,r1vers_2) + ...
                 cos(beta2)*hvers_2);

%vinf @ Mars = 0;
             
[ m2, T2 ] = Conway(TOF2, N_rev2, q2, r1norm_2, r2norm_2, r1vers_2, r2vers_2, hvers_2, hh_2, v1_2, v2_2, muS, data_stacks);

%thrid launch
[kepEarth_3, ~] = uplanet(t03,3);
[kepMars_3, ~]    = uplanet(t03 + TOF3,4);

[R1_3, v1_3] = kep2car2(kepEarth_3, muS); %km....
[R2_3, v2_3] = kep2car2(kepMars_3, muS);  %km....

%definition of plane of motion
r1norm_3 = norm(R1_3);
r2norm_3 = norm(R2_3);

r1vers_3 = R1_3/r1norm_3;
r2vers_3 = R2_3/r2norm_3;

hh_3 = cross(r1vers_3, r2vers_3);
hvers_3 = hh_3/norm(hh_3);

if hh_3(3) <0
    hvers_3 = -hvers_3;
end

%adding TMI maneuver
v1_3 = v1_3 + v_inf3*(sin(beta3)*cos(alpha3)*r1vers_3 + ...
                 sin(beta3)*sin(alpha3)*cross(hvers_3,r1vers_3) + ...
                 cos(beta3)*hvers_3);

%vinf @ Mars = 0;
             
[ m3, T3 ] = Conway(TOF3, N_rev3, q3, r1norm_3, r2norm_3, r1vers_3, r2vers_3, hvers_3, hh_3, v1_3, v2_3, muS, data_stacks);

if ~all(isnan([m1, m2, m3]))
    t1 = max(abs(T1)); t2 = max(abs(T2)); t3 = max(abs(T3));
    y = t1 + t2 + t3;
else 
    y = inf;
end

end

