
function [r, TH, L, gamma1, gamma2, gamma, a_inplane, vr, vt, v1tra, v2tra, vnorm, time, dmdt, m, T_inplane, TOFr] = Conway(TOF, N_rev, data_stacks)
[kepEarth, muS] = uplanet(data_stacks.t0sym, 3);
data_stacks.tmax = data_stacks.t0sym + TOF;
[kepMars, ~] = uplanet(data_stacks.tmax,4);
[r1, v1] = kep2car2(kepEarth, muS); %km....
[r2, v2] = kep2car2(kepMars, muS);  %km....

%definition of plane of motion
r1norm = norm(r1);
r2norm = norm(r2);

r1vers =r1/r1norm;
r2vers =r2/r2norm;

RIvcRFv = cross(r1vers, r2vers);
RCRRv = RIvcRFv/norm(RIvcRFv);
rivDrfv = dot(r1vers, r2vers);

if RIvcRFv(3) <0
    RCRRv = -RCRRv;
end

%decomposition of velocity
v1perp = dot(v1,RCRRv);
v1plan = v1 - v1perp*RCRRv;
v2perp = dot(v2,RCRRv);
v2plan = v2 - v2perp*RCRRv;

if norm(RIvcRFv)<1e-10
    if rivDrfv > 0
        phi = 0;
    else 
        phi = pi;
    end
else
    if  RIvcRFv(3)> 0
        phi = acos(rivDrfv);
    else 
        phi = 2*pi -  acos(rivDrfv);
    end
end

L = phi + 2*pi*N_rev;

v1rad = r1vers*dot(v1plan,r1vers);  v1tra = v1plan - v1rad;
v1rad = r2vers*dot(v2plan,r2vers);  v2tra = v2plan - v1rad;

gamma1 = atan(dot(v1plan, r1vers)/norm(v1tra));
gamma2 = atan(dot(v2plan, r2vers)/norm(v2tra));

th_dot1 = norm(v1tra)/r1norm;
th_dot2 = norm(v2tra)/r2norm;

%Conway algorithm 
a = 1/r1norm;
b = -tan(gamma1)/r1norm;
c = 1/(2*r1norm) * ( (muS/(r1norm^3 * th_dot1^2)) -1 );
ABC = [a b c];

opts = optimset('TolX', 1e-14, 'TolFun', 1e-14, 'MaxFunEvals', 1e3);
[dd, feval] = fzero(@(d) find_d(d, ABC, muS, L, r2, th_dot2, gamma2, TOF), 0, opts);

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
par.muS = muS;

A0 = [30*L^2   -10*L^3   L^4; ...
     -48*L      18*L^2  -2*L^3; ...
      20       -8*L      L^2];

B0 = [1/r2norm - (a + b*L + c*L^2 + dd*L^3); ...
    -tan(gamma2)/r2norm-(b +2*c*L + 3*dd*L^2); ...
     muS/(r2norm^4 * th_dot2^2) - (1/r2norm + 2*c + 6*dd*L)];

e = 0.5/L^6*A0(1,:)*B0;
f = 0.5/L^6*A0(2,:)*B0;
g = 0.5/L^6*A0(3,:)*B0;

% % Out of plane
% a_out = 0;
% b_out = v1(3)/th_dot1;
% q = 0.1;
% A1 = [q*L L^2; -(q-1) L];
% B1 = [r2(3) - a_out - b_out*L; ...
%       (v2(3)/th_dot2)- b_out];
% CD = 1/(L^q)*A1*B1;
% c_out = CD(1); d_out = CD(2);

[~, ~, TH] = TOFcomp(dd, ABC, muS, L, r2, th_dot2, gamma2);

%computation of vector of time of flight
n_int = 10000; d = dd;
for i = 1:n_int
    tm(i) = 0.5*(TH(i+1)+TH(i));
end

r = 1./(a + b*TH + c*TH.^2 + d*TH.^3 + e*TH.^4 + f*TH.^5 + g*TH.^6);
rm = 1./(a + b*tm + c*tm.^2 + d*tm.^3 + e*tm.^4 + f*tm.^5 + g*tm.^6);
dTOF = (abs(r.^4/muS .* (1./r + 2*c + 6*d*TH + 12*e*TH.^2 + 20*f*TH.^3 + 30*g*TH.^4)) ).^0.5;
dTOFm = (abs(rm.^4/muS .* (1./rm + 2*c + 6*d*tm + 12*e*tm.^2 + 20*f*tm.^3 + 30*g*tm.^4)) ).^0.5;

delta_TH = TH(2) - TH(1);

%integration
TOFr = zeros(n_int + 1,1);
for i = 2:n_int+1
    TOFr(i) = TOFr(i-1) + dTOF(i-1) + 4*dTOFm(i-1) + dTOF(i); % [s]
end
TOFr = TOFr*delta_TH/6/86400;                                 % [days] 

theta_dot = ((muS./r.^4)./(1./r + 2*c + 6*d*TH + 12*e*TH.^2 +20*f*TH.^3 + 30*g*TH.^4)).^0.5;

if isreal(theta_dot)
    
    %flight path angle
    gamma = atan(-r.*(b + 2*c*TH + 3*d*TH.^2 + 4*e*TH.^3 + 5*f*TH.^4 + 6*g*TH.^5));
    
    
    a_inplane = -(muS./(2*r.^3.*cos(gamma))) .* (6*d + 24*e*TH + 60*f*TH.^2 + 120*g*TH.^3 - tan(gamma)./r)./ ...
                  (1./r + 2*c + 6*d*TH + 12*e*TH.^2 + 20*f*TH.^3 + 30*g*TH.^4).^2;
              
    %velocity components
    vr = -r.^2.*theta_dot.*(b + 2*c*TH +3*d*TH.^2 + 4*e*TH.^3 + 5*f*TH.^4 + 6*g*TH.^5);
    vt = r.*theta_dot;
    vnorm = sqrt(vt.^2 + vr.^2);
    
    %propellant mass and thrust required
    time = TOFr*86400;
    m = zeros(size(TH)); m(end) = data_stacks.Mdry;
       
    %integration
    for i = length(m):-1:2

        hhh = time(i-1) - time(i);
        K1 = massflowrate(time(i)       , m(i)            , a_inplane(i)                       , data_stacks);
        K2 = massflowrate(time(i)+ hhh/2, m(i) + hhh/2*K1 , 0.5*(a_inplane(i) + a_inplane(i-1)), data_stacks);
        K3 = massflowrate(time(i)+ hhh/2, m(i) + hhh/2*K2 , 0.5*(a_inplane(i) + a_inplane(i-1)), data_stacks);
        K4 = massflowrate(time(i)+ hhh  , m(i) + hhh  *K3 , a_inplane(i-1)                     , data_stacks);
        
        m(i-1) = m(i) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end
    dmdt = massflowrate(time, m, a_inplane, data_stacks);
    T_inplane = a_inplane.*m *1000;
else
    r = NaN;
    TH = NaN;
    L = NaN; 
    gamma1 = NaN;
    gamma2 = NaN;
    v1tra = NaN;
    v2tra = NaN; 
    vnorm = NaN;
    TOFr = NaN;
    gamma = NaN;
    a_inplane = NaN;
    vr = NaN;
    vt = NaN;
    time = NaN;
    dmdt = NaN;
    m = NaN;
    T_inplane = NaN;
end
% gamma, a_inplane, vr, vt, time, dmdt, m, T_inplane
end