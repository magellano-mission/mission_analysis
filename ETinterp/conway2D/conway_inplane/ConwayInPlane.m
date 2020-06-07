function [T, m, r, TH] = ConwayInPlane(TOF, N_rev, r1norm, r2norm, r1vers, r2vers, hvers, v1, v2, miS, data_stacks)

n_int = data_stacks.n_int;
%decomposition of velocity

phi = acos( max( min( dot(r1vers, r2vers), 1) , -1) );

if  hvers(3) <= 0
    phi = 2*pi - phi;
end

L = phi + 2*pi*N_rev;

v1r = r1vers*dot(v1, r1vers);  v1t = v1 - v1r;
v2r = r2vers*dot(v2, r2vers);  v2t = v2 - v2r;

gamma1 = atan2(norm(v1r), norm(v1t));
gamma2 = atan2(norm(v2r), norm(v2t));

th_dot1 = norm(v1t)/r1norm;
th_dot2 = norm(v2t)/r2norm;

%Conway algorithm
a = 1 / r1norm;
b = -tan(gamma1) / r1norm;
c = 1 / (2*r1norm) * ( (miS/ (r1norm^3 * th_dot1^2)) -1 );
ABC = [a b c];

% opts = optimset('Display','off');
[d_plan, ~, ~] = fzero(@(d) find_d(d, ABC, miS, L, r2norm, th_dot2, gamma2, TOF, n_int), 1e-10);
% flag;

A0 = [30*L^2   -10*L^3   L^4; ...
    -48*L      18*L^2  -2*L^3; ...
    20       -8*L      L^2];

B0 = [1/r2norm - (a + b*L + c*L^2 + d_plan*L^3); ...
    -tan(gamma2)/r2norm - (b +2*c*L + 3*d_plan*L^2); ...
    miS/(r2norm^4 * th_dot2^2) - (1/r2norm + 2*c + 6*d_plan*L)];

EFG = 0.5/L^6 * A0 * B0;
e = EFG(1);
f = EFG(2);
g = EFG(3);

[~, r, TH] = TOFcomp(d_plan, ABC, miS, L, r2norm, th_dot2, gamma2, n_int);

%TH: [n_int,1]
TH2 = TH.^2;
TH3 = TH.^3;
TH4 = TH.^4;
TH5 = TH.^5;

% (miS./r.^4)./(1./r + 2*c + 6*d_plan*TH + 12*e*TH2 +20*f*TH3 + 30*g*TH4);
theta_dot = ((miS./r.^4)./(1./r + 2*c + 6*d_plan*TH + 12*e*TH2 +20*f*TH3 + 30*g*TH4)).^0.5;

if isreal(theta_dot)
    
    %flight path angle
    gamma = atan(-r.*(b + 2*c*TH + 3*d_plan*TH2 + 4*e*TH3 + 5*f*TH4 + 6*g*TH5));
    
    acc = -(miS./(2*r.^3.*cos(gamma))) .* (6*d_plan + 24*e*TH + 60*f*TH2 + 120*g*TH3 - tan(gamma)./r)./ ...
        (1./r + 2*c + 6*d_plan*TH + 12*e*TH2 + 20*f*TH3 + 30*g*TH4).^2;
    
    d_time = 1./theta_dot;
    dl = TH(2) - TH(1);
    time = zeros(1,n_int);
    for i = 2:n_int-1
        time(i) = time(i-1) + dl/6 * (d_time(i-1) + 4*d_time(i) + d_time(i + 1));
    end
    time(end) = time(end-1) + d_time(end-1)*dl;
    
    %propellant mass and thrust required
    m = zeros(size(TH));        m(end) = data_stacks.Mdry; %[kg]
    %RK4 backwards integration
    for i = length(m):-1:2
        
        hhh = time(i-1) - time(i);
        K1 = massflowrate(time(i)       , m(i)            , acc(i)                 ,     data_stacks);
        K2 = massflowrate(time(i)+ hhh/2, m(i) + hhh/2*K1 , 0.5*(acc(i) + acc(i-1)),     data_stacks);
        K3 = massflowrate(time(i)+ hhh/2, m(i) + hhh/2*K2 , 0.5*(acc(i) + acc(i-1)),     data_stacks);
        K4 = massflowrate(time(i)+ hhh  , m(i) + hhh  *K3 , acc(i-1)               ,     data_stacks);
        
        m(i-1) = m(i) + hhh/6*(K1 + 2*K2 + 2*K3 + K4);
    end
    
    T = acc.*m*1e3;
else
    m = inf;
    T = inf;
end
