function [y] = conway_flyby(x,data_stacks)
t0              = x(1); 
TOF12           = x(2); 
TOF23           = x(3);
N_rev1          = round(x(4));  
N_rev2          = round(x(5));
q1              = x(6);     
q2              = x(7);
%departure
v_inf           = x(8);
alpha           = x(9); 
beta            = x(10);
%gravity assist
v_infga         = x(11);
alphagam        = x(12); 
alphagap        = x(13); 
betagam         = x(14); 
betagap         = x(15); %flyby

[kep2, muS]   = uplanet(t0 + TOF12, 3); %flyby at Earth
[R2, v2] = kep2car2(kep2, muS);  %km....

muE = 3.986004330000000e+05;
rE = 6.371010000000000e+03;

%gravity assist check
xversga = v2/norm(v2);
zversga = cross(R2, v2)/norm(cross(R2,v2));
yversga = cross(xversga, zversga);

V_infm = (cos(betagam)* (cos(alphagam)*xversga + sin(alphagam)*yversga) + sin(betagam)*zversga);
V_infp = (cos(betagap)* (cos(alphagap)*xversga + sin(alphagap)*yversga) + sin(betagap)*zversga);
delta = acos(dot(V_infm, V_infp));
e = 1/ sin(delta/2);
rp = muE/norm(v_infga)^2 *(e-1);

if rp < 1.3*rE
        y = 999999*(1*rand(2,1));
else
    
    delta_vga = v_infga*(V_infp - V_infm);
    v2 = v2 + delta_vga;

    [kep1, ~] = uplanet(t0              ,   3);
    [kep3, ~]   = uplanet(t0 + TOF12 + TOF23,   4);

    [R1, v1] = kep2car2(kep1, muS); %km....
    [R3, v3] = kep2car2(kep3, muS);  %km....

    %definition of plane of motion
    r1norm = norm(R1);
    r2norm = norm(R2);
    r3norm = norm(R3);

    r1vers = R1/r1norm;
    r2vers = R2/r2norm;
    r3vers = R3/r3norm;

    RIvcRFv12 = cross(r1vers, r2vers);
    RCRRv12 = RIvcRFv12/norm(RIvcRFv12);

    if RIvcRFv12(3) <0
        RCRRv12 = -RCRRv12;
    end

    RIvcRFv23 = cross(r2vers, r3vers);
    RCRRv23 = RIvcRFv23/norm(RIvcRFv23);

    if RIvcRFv23(3) <0
        RCRRv23 = -RCRRv23;
    end

    [ m2, T2 ] = Conway(TOF23, N_rev2, q2, r2norm, r3norm, r2vers, r3vers, RCRRv23, RIvcRFv23, v2, v3, muS, data_stacks);

    data_stacks.Mdry = m2(end);

    %TMI maneuver
    V_infdep = v_inf*(sin(beta)*cos(alpha)*r1vers + ...
                     sin(beta)*sin(alpha)*cross(RCRRv12,r1vers) + ...
                     cos(beta)*RCRRv12);
    v1 = v1 + V_infdep;

    [ m1, T1 ] = Conway(TOF12, N_rev1, q1, r1norm, r2norm, r1vers, r2vers, RCRRv12, RIvcRFv12, v1, v2, muS, data_stacks);

    %computation of fly-by

    if ~isnan(m1)
        if ~isnan(m2)
            y(1) = max(max(abs(T1)),max(abs(T2)));
            y(2) = m1(1) - m2(end);
        else
            y = 999999*(1*rand(2,1));
        end
    else
            y = 999999*(1*rand(2,1));
    end
    end
end

