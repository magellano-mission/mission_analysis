%% porkchop plot ET
clear all, close all, clc
%simulation parameters
data_stacks.Isp = 4300;                                      % specific impulse [s]
data_stacks.Mdry = 1600;                                      % Total Mass of the s/c [kg]
data_stacks.n_int = 1000;
v_inf      =  1.019932213771344;
alpha      =  1.339209562214939;
beta       =          2.689364771171173;
% t0 = 9.618881860211146e+03; 
TOFrange = (700:5:1300);
t0 = date2mjd2000([2024 1 1 0 0 0]); 
tf = date2mjd2000([2029 1 1 0 0 0]); 
t0range = (t0:5:tf);
N_rev = 1;  q = 5.14;
  
% v_inf = 1.026859552871398; alpha = 0.824387529700576; beta = 0.450043526794778;

M = zeros(length(TOFrange), length(t0range)); T = M;
for j = 1:length(t0range)
for i = 1:length(TOFrange)
    wbb = waitbar((j)/(length(t0range)));
    t0 = t0range(j);
    TOF = TOFrange(i);

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

    [ m, t ] = Conway(TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, RCRRv, RIvcRFv, v1, v2, muS, data_stacks);
    
        if isnan(m) 
        M(i,j) = NaN; 
        T(i,j) = NaN;
        else
        M(i,j) = m(1) - m(end); 
        T(i,j) = max(abs(t));
        end

end
end

delete(wbb)

%% porkchop
open('porkchop1.mat')
% M = ans.M;
% T = ans.T;
% t0range = ans.t0range;
% TOFrange = ans.TOFrange;
Mfiltered = M;
Tfiltered = T;

%processing
for j = 1:length(t0range)
for i = 1:length(TOFrange)
    
    if T(i,j)>0.2
        Mfiltered(i,j) = NaN; 
        Tfiltered(i,j) = NaN;
    end
end
end
figure() 
subplot(2,1,1)

s1 = pcolor( t0range, TOFrange, Mfiltered ); hold on
s1.FaceColor = 'Interp';
s1.EdgeColor = 'Interp'; grid on
xlabel('Departure Date [MJD2000]', 'Interpreter','Latex'), ylabel('TOF [days]','Interpreter','Latex')
title('Required propellant mass','Interpreter','Latex')
colorbar
subplot(2,1,2)

s2 = pcolor(t0range, TOFrange, Tfiltered ); hold on
s2.FaceColor = 'Interp';
s2.EdgeColor = 'Interp'; grid on
title('Maximum Thrust Required','Interpreter','Latex')
xlabel('Departure Date [MJD2000]','Interpreter','Latex'), ylabel('TOF [days]','Interpreter','Latex')
colorbar