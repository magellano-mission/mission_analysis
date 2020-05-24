%% porkchop plot ET
%simulation parameters
data_stacks.Isp = 5000;                                      % specific impulse [s]
data_stacks.Mdry = 2800;                                      % Total Mass of the s/c [kg]
data_stacks.n_int = 1000;

TOFrange = (700:10:1500);
t0 = mjd2date([2024 1 1 0 0]); 
tf = mjd2date([2029 1 1 0 0]); 
t0range = (t0:10:tf);
N_rev = 2;  q = 6.28;
  
v_inf = 0; alpha = 0; beta = 0;
v_infcap = 0; alphacap = 0; betacap = 0;

M = zeros(length(TOFrange), length(t0range)); T = M;
for j = 1:length(t0range)
for i = 1:length(TOFrange)
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

    %adding v_inf at Mars capture
    v2 = v2 + v_infcap*(sin(betacap)*cos(alphacap)*r2vers + ...
                     sin(betacap)*sin(alphacap)*cross(RCRRv,r2vers) + ...
                     cos(betacap)*RCRRv);

    [ m, T ] = Conway(TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, RCRRv, RIvcRFv, v1, v2, muS, data_stacks);
    
    
    if ~isnan(m) 
        M(i,j) = NaN; 
        T(i,j) = NaN;
    elseif max(abs(T))<0.25
        M(i,j) = NaN; 
        T(i,j) = NaN;
    else
        M(i,j) = m(1) - m(end); 
        T(i,j) = max(abs(T));
    end
end
end
%porkchop
figure() 
sgtitle(' Required propellant mass')
s1 = pcolor(t0range, TOFrange, M ); hold on
s1.FaceColor = 'Interp';
s1.EdgeColor = 'Interp';

figure()
sgtitle(' Required trhrust')
s2 = pcolor(t0range, TOFrange, T ); hold on
s2.FaceColor = 'Interp';
s2.EdgeColor = 'Interp';