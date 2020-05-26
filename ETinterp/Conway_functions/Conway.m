function [ m, T, r, z, s, vr, vt, vz, acc_inplane, acc_out, acc, TH, L, gamma1, gamma2, gamma, v1perp, v2perp, v1tra, v2tra, vnorm, dmdt, T_inplane, T_outplane, theta_dot, time, TOFr] = ...
    Conway(TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, hvers, hh, v1, v2, muS, data_stacks)

n_int = data_stacks.n_int;            
%decomposition of velocity
v1perp = dot(v1,hvers);
v2perp = dot(v2,hvers);

v1plan = v1 - v1perp*hvers;
v2plan = v2 - v2perp*hvers;

phi = acos( max( min( dot(r1vers, r2vers), 1) , -1) );

if  hh(3) <= 0
    phi = 2*pi - phi;
end

L = phi + 2*pi*N_rev;

v1rad = r1vers*dot(v1plan,r1vers);  v1tra = v1plan - v1rad;
v2rad = r2vers*dot(v2plan,r2vers);  v2tra = v2plan - v2rad;

gamma1 = atan( dot(v1plan, r1vers)/norm(v1tra) );
gamma2 = atan( dot(v2plan, r2vers)/norm(v2tra) );

th_dot1 = norm(v1tra)/r1norm;
th_dot2 = norm(v2tra)/r2norm;

%Conway algorithm 
a_plan = 1 / r1norm;
b_plan = -tan(gamma1) / r1norm;
c_plan = 1 / (2*r1norm) * ( (muS/ (r1norm^3 * th_dot1^2)) -1 );
ABC = [a_plan b_plan c_plan];

fun = @(d) find_d(d, ABC, muS, L, r2norm, th_dot2, gamma2, TOF, n_int);
opts = optimset('Display','off');
d_plan = fzero(fun, 1e-10, opts);

A0 = [30*L^2   -10*L^3   L^4; ...
     -48*L      18*L^2  -2*L^3; ...
      20       -8*L      L^2];

B0 = [1/r2norm - (a_plan + b_plan*L + c_plan*L^2 + d_plan*L^3); ...
     -tan(gamma2)/r2norm - (b_plan +2*c_plan*L + 3*d_plan*L^2); ...
     muS/(r2norm^4 * th_dot2^2) - (1/r2norm + 2*c_plan + 6*d_plan*L)];

EFG = 0.5/L^6 * A0 * B0;
e_plan = EFG(1); f_plan = EFG(2); g_plan = EFG(3);

[~, r, TH] = TOFcomp(d_plan, ABC, muS, L, r2norm, th_dot2, gamma2, n_int);

%TH: [n_int,1]
TH2 = TH.^2; TH3 = TH.^3;
TH4 = TH.^4; TH5 = TH.^5;

theta_dot = ((muS./r.^4)./(1./r + 2*c_plan + 6*d_plan*TH + 12*e_plan*TH2 +20*f_plan*TH3 + 30*g_plan*TH4)).^0.5;

if isreal(theta_dot)   
    
    %flight path angle
    gamma = atan(-r.*(b_plan + 2*c_plan*TH + 3*d_plan*TH2 + 4*e_plan*TH3 + 5*f_plan*TH4 + 6*g_plan*TH5));
    
    acc_inplane = -(muS./(2*r.^3.*cos(gamma))) .* (6*d_plan + 24*e_plan*TH + 60*f_plan*TH2 + 120*g_plan*TH3 - tan(gamma)./r)./ ...
                  (1./r + 2*c_plan + 6*d_plan*TH + 12*e_plan*TH2 + 20*f_plan*TH3 + 30*g_plan*TH4).^2;
              
    %velocity components
    vr = -r.^2.*theta_dot.*(b_plan + 2*c_plan*TH +3*d_plan*TH2 + 4*e_plan*TH3 + 5*f_plan*TH4 + 6*g_plan*TH5);
    vt = r.*theta_dot;
    
    TDD1 = 4*tan(gamma)./(1./r + 2*c_plan + 6*d_plan*TH + 12*e_plan*TH2 + 20*f_plan*TH3 + 30*g_plan*TH4);
    TDD2 = (6*d_plan + 24*e_plan*TH + 60*f_plan*TH2 + 120*g_plan*TH3 - tan(gamma)./r)/ ...
           (1./r + 2*c_plan + 6*d_plan*TH + 12*e_plan*TH2 + 20*f_plan*TH3 + 30*g_plan*TH4).^2;
    
    th_dot_dot = -muS./(2*r.^4).*(TDD1+TDD2);
    d_time = 1./theta_dot;
    dl = TH(2) - TH(1);
    time = zeros(1,n_int);
    for i = 2:n_int-1
        time(i) = time(i-1) + dl/6 * (d_time(i-1) + 4*d_time(i) + d_time(i + 1));
    end
    time(end) = time(end-1) + d_time(end-1)*dl;

%     dTOF  = (abs(r.^4/muS .* (1./r + 2*c_plan + 6*d_plan*TH + 12*e_plan*TH2 + 20*f_plan*TH3 + 30*g_plan*TH4)) ).^0.5;
% 
%     hh = TH(2) - TH(1);
%     %Cavalieri method
%     time = zeros(size(TH));
%     for i = 2:n_int  
%         tm = 0.5*(TH(i-1)+TH(i));   
%         rm = 1./(a_plan + b_plan*tm + c_plan*tm^2 + d_plan*tm^3 + e_plan*tm^4 + f_plan*tm^5 + g_plan*tm^6);   
%         dTOFm = ( abs(rm^4/muS * (1/rm + 2*c_plan + 6*d_plan*tm + 12*e_plan*tm^2 + 20*f_plan*tm^3 + 30*g_plan*tm^4)) )^0.5;
%         time(i) = time(i-1) + hh/6*(dTOF(i-1) + 4*dTOFm + dTOF(i));
%     end
    TOFr = time/86400; %[days]
    
    z1 = 0;    z2 = 0;
    zd1 = v1perp;    zd2 = v2perp;
    
    % Out of plane
    a_out = z1;
    b_out = zd1/th_dot1;
    
    A1 = [q*L -L^2; -(q-1) L];
    
    B1 = [z2 - a_out - b_out*L; ...
         (zd2/th_dot2)- b_out]; %?

    CD = 1/(L^q)*A1*B1;
    c_out = CD(1); d_out = CD(2);
    
    z = a_out + b_out*TH + c_out*TH.^(q-1) + d_out*TH.^q; %[km]
    vz = (b_out + c_out*(q-1)*TH.^(q-2) + d_out*q*TH.^(q-1)).*theta_dot;
    vnorm = (vt.^2 + vr.^2 + vz.^2).^0.5;
    
    s = (r.^2 + z.^2).^0.5;
    
    acc_out = muS./s.^3 .*z + (c_out*(q-1)*(q-2)*TH.^(q-3) + d_out*q*(q-1)*TH.^(q-2)).*theta_dot.^2 + ...
                (b_out + c_out*(q-1)*TH.^(q-2) + d_out*q*TH.^(q-1)).*th_dot_dot; %[km/s^2]
    
    %propellant mass and thrust required
    m = zeros(size(TH));        m(end) = data_stacks.Mdry; %[kg]
    acc = (acc_inplane.^2 + acc_out.^2).^0.5;  
    %RK4 backwards integration
    for i = length(m):-1:2

        hhh = time(i-1) - time(i);
        K1 = massflowrate(time(i)       , m(i)            , acc(i)                 ,     data_stacks);
        K2 = massflowrate(time(i)+ hhh/2, m(i) + hhh/2*K1 , 0.5*(acc(i) + acc(i-1)),     data_stacks);
        K3 = massflowrate(time(i)+ hhh/2, m(i) + hhh/2*K2 , 0.5*(acc(i) + acc(i-1)),     data_stacks);
        K4 = massflowrate(time(i)+ hhh  , m(i) + hhh  *K3 , acc(i-1)               ,     data_stacks);
        
        m(i-1) = m(i) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end
    dmdt = massflowrate(time, m, acc, data_stacks);
    T_inplane = acc_inplane .*m * 1000;
    T_outplane = acc_out .*m * 1000;
    T = (T_inplane.^2 + T_outplane.^2).^0.5;
    else
    r = NaN; s=r; z =r;
    TH = NaN;
    L = NaN; 
    gamma1 = NaN;
    gamma2 = NaN;
    gamma = NaN;
    v1perp = NaN;
    v2perp = NaN;
    v1tra = NaN;
    v2tra = NaN; 
    vnorm = NaN;
    TOFr = NaN;
    acc = NaN; acc_inplane = NaN; acc_out = NaN;
    vr = NaN;
    vt = NaN;
    vz = NaN;
    time = NaN;
    dmdt = NaN;
    m = NaN;
    T = NaN; T_inplane = NaN; T_outplane = NaN;
end
% gamma, a_inplane, vr, vt, time, dmdt, m, T_inplane
end