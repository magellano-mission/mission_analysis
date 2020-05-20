%Conway-based approach

%definition of parameters
close all, clear, clc
% Figure Initialization    
load('MagellanoColorMap.mat');
DefaultOrderColor = get(0, 'DefaultAxesColorOrder');
NewOrderColor = [0.9490    0.4745    0.3137
                 0.1020    0.6667    0.74120
                 155/255   155/255   155/255
                 DefaultOrderColor];  
             
set(0,'DefaultFigureColormap', MagellanoColorMap);
set(0, 'DefaultAxesColorOrder', NewOrderColor);
set(0,'DefaultLineLineWidth', 2)
set(0,'DefaultLineMarkerSize', 10)
set(0, 'DefaultFigureUnits', 'normalized');
set(0, 'DefaultFigurePosition', [0 0 1 1]);
set(0, 'DefaultTextFontSize', 18);
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultAxesXGrid', 'on')
set(0, 'DefaultAxesYGrid', 'on')
set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');

data_stacks.t0sym = date2mjd2000([2024 1 1 0 0 0]);
data_stacks.tmax = date2mjd2000([2030 1 1 0 0 0]);

data_stacks.isPerturbed = 0;
data_stacks.isInterp = 1;
data_stacks.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
data_stacks.Isp = 3000;                                      % specific impulse [s]
data_stacks.M0 = 500;                                      % Total Mass of the s/c [kg]
data_stacks.Across_sun = 10;                                % Cross area related to the sun [m^2]

data_stacks.event = 1;
data_stacks.opt = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'InitialStep', 1e-12, 'Event', @event_interp_LT);

Thrust = 0.2;                                               % value of thrust [N]
data_stacks.T = [Thrust; 0; 0];                           % thrust [N] (@TNH) (constant profile for the moment)

%definition of interplanetary arc
TOF = 5*365;%[days]
N_rev= 3;
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

vErad = r1vers*dot(v1plan,r1vers);  v1tra = v1plan - vErad;
vMrad = r2vers*dot(v2plan,r2vers);  v2tra = v2plan - vMrad;

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

[I, ~, TH] = TOFcomp(dd, ABC, muS, L, r2, th_dot2, gamma2);

n_int = 10000; d = dd;
for i = 1:n_int
    tm(i) = 0.5*(TH(i+1)+TH(i));
end
r = 1./(a + b*TH + c*TH.^2 + d*TH.^3 + e*TH.^4 + f*TH.^5 + g*TH.^6);
rm = 1./(a + b*tm + c*tm.^2 + d*tm.^3 + e*tm.^4 + f*tm.^5 + g*tm.^6);
dTOF = (abs(r.^4/muS .* (1./r + 2*c + 6*d*TH + 12*e*TH.^2 + 20*f*TH.^3 + 30*g*TH.^4)) ).^0.5;
dTOFm = (abs(rm.^4/muS .* (1./rm + 2*c + 6*d*tm + 12*e*tm.^2 + 20*f*tm.^3 + 30*g*tm.^4)) ).^0.5;
hh = TH(2) - TH(1);
TOFr = zeros(n_int + 1,1);
for i = 2:n_int+1
    TOFr(i) = TOFr(i-1) + dTOF(i-1) + 4*dTOFm(i-1) + dTOF(i); %s
end
TOFr = TOFr*hh/6/86400; %days 
RE = zeros(length(TOFr), 3); RM = RE;
REnorm = zeros(length(TOFr),1); RMnorm = REnorm;
for i =1:length(TOFr)
    kepE = uplanet(data_stacks.t0sym+TOFr(i), 3);
    kepM = uplanet(data_stacks.t0sym+TOFr(i), 4);
    RE(i,:) = kep2car2(kepE, muS);
    RM(i,:) = kep2car2(kepM, muS);
    REnorm(i) = norm(RE(i,:));
    RMnorm(i) = norm(RM(i,:));
end
figure()
sgtitle('Sun distance')
plot(TOFr,RMnorm, 'DisplayName', 'Mars'), hold on
plot(TOFr,REnorm, 'DisplayName', 'Earth'), hold on
plot(TOFr, r, 'DisplayName', 's/c'), hold on, legend()

ABCDEFG = [a b c d e f g];

theta_dot = ((muS./r.^4)./(1./r + 2*c + 6*d*TH + 12*e*TH.^2 +20*f*TH.^3 + 30*g*TH.^4)).^0.5;

if isreal(theta_dot)
    gamma = atan(-r.*(b + 2*c*TH + 3*d*TH.^2 + 4*e*TH.^3 + 5*f*TH.^4));
    a_inplane = -muS./(2*r.^3.*cos(gamma)) .* (6*d + 24*e*TH + 60*f*TH.^2 + 120*g*TH.^3 - tan(gamma)./r)./ ...
                  (1./r + 2*c + 6*d*TH + 12*e*TH.^2 + 20*f*TH.^3 + 30*g*TH.^4).^2;
    vr = -r.^2.*theta_dot.*(b + 2*c*TH +3*d*TH.^2 + 4*e*TH.^3 + 5*f*TH.^4 + 6*g*TH.^5);

    time = TOFr*86400;
    m = zeros(size(TH)); m(end) = 1000;
    %integration
    for i = length(m):-1:2
        
        hhh = time(i-1) - time(i);
        K1 = massflowrate(time(i)       , m(i)          , a_inplane(i)                       , data_stacks);
        K2 = massflowrate(time(i)+ hhh/2, m(i)+hhh/2*K1 , 0.5*(a_inplane(i) + a_inplane(i-1)), data_stacks);
        K3 = massflowrate(time(i)+ hhh/2, m(i)+hhh/2*K2 , 0.5*(a_inplane(i) + a_inplane(i-1)), data_stacks);
        K4 = massflowrate(time(i-1)     , m(i)+hhh*K3   , a_inplane(i-1)                     , data_stacks);
        
        m(i-1) = m(i) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end
    dmdt = massflowrate(time, m, a_inplane, data_stacks);
    T_inplane = a_inplane.*m *1000;
else
    gamma = NaN;
    a_inplane = NaN;
    vr = NaN;
    dmdt = NaN;
end

figure()
sgtitle('Thrust Profile')
subplot(3,1,1), plot(time, a_inplane), title('a_{inplane}')
subplot(3,1,2), plot(time, T_inplane), title('inplane T')
subplot(3,1,3), plot(time, m), title('mass')

%% functions

function [dmdt] = massflowrate(time ,m, a, data_stacks)
T = a.*m*1000;
dmdt = -abs(T)/(9.81*data_stacks.Isp);
end
function errorTOF = find_d(d, ABC, muS, th_f, rM, th_dotM, gammaM, TOF)
I = TOFcomp(d, ABC, muS, th_f, rM, th_dotM, gammaM);
errorTOF = TOF*86400 - I;
end
function [I, r, TH] = TOFcomp(d, ABC, muS, th_f, rM, th_dotM, gammaM)
a = ABC(1); b = ABC(2); c = ABC(3);
rMnorm = norm(rM);
A0 = [30*th_f^2   -10*th_f^3   th_f^4; ...
     -48*th_f    18*th_f^2  -2*th_f^3; ...
      20        -8*th_f      th_f^2];

A1 = [1/rMnorm-(a + b*th_f + c*th_f^2 + d*th_f^3); ...
    -tan(gammaM)/rMnorm-(b +2*c*th_f + 3*d*th_f^2); ...
     muS/(rMnorm^4*th_dotM^2)-(1/rMnorm + 2*c + 6*d*th_f)];
   
e = 0.5/th_f^6*A0(1,:)*A1;
f = 0.5/th_f^6*A0(2,:)*A1;
g = 0.5/th_f^6*A0(3,:)*A1;

n_int = 10000;
TH = linspace(0, th_f, n_int + 1);
for i = 1:n_int
    tm(i) = 0.5*(TH(i+1)+TH(i));
end

r = 1./(a + b*TH + c*TH.^2 + d*TH.^3 + e*TH.^4 + f*TH.^5 + g*TH.^6);
rm = 1./(a + b*tm + c*tm.^2 + d*tm.^3 + e*tm.^4 + f*tm.^5 + g*tm.^6);

dTOF = (abs(r.^4/muS .* (1./r + 2*c + 6*d*TH + 12*e*TH.^2 + 20*f*TH.^3 + 30*g*TH.^4)) ).^0.5;
dTOFm = (abs(rm.^4/muS .* (1./rm + 2*c + 6*d*tm + 12*e*tm.^2 + 20*f*tm.^3 + 30*g*tm.^4)) ).^0.5;

hh = TH(2) - TH(1);
%Cavalieri method
I = 0;
for i = 2:n_int+1
    I = I + dTOF(i-1) + 4*dTOFm(i-1) + dTOF(i);
end
I = I*hh/6;
end