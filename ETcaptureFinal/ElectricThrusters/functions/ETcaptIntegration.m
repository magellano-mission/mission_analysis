function dY = ETcaptIntegration(t, Y, data)
%{
ODE function to integrate the Thrusted electric arches
%}

% States
R = Y(1:3);
V = Y(4:6);
M = Y(7);

% input checks
% if isrow(R)
%     R = R';
% end
% if isrow(V)
%     V = V';
% end
%% Gravity perturb
R_pl = data.R_pl;                          % planet radius
J2 = data.J2;                              % J2 perturbation constant
J3 = data.J3;                              % J2 perturbation constant
J4 = data.J4;                              % J2 perturbation constant
mi = data.mi;                              % earth gravitational constant
w = data.wM;                                
r = norm(R);


% Rotation matrix from ECI to ECEF frame
B = [cos(w*t) sin(w*t) 0 ;
    -sin(w*t) cos(w*t) 0 ;
        0         0      1];

% S/c position in ECEF frame
s_sc = B*R;

% J2 Effect
k = 3/2 * J2 * mi * R_pl^2 / r^4;

a_j2_x = k*(s_sc(1)/r *(5 *s_sc(3)^2/r^2 -1));
a_j2_y = k*(s_sc(2)/r *(5 *s_sc(3)^2/r^2 -1));
a_j2_z = k*(s_sc(3)/r *(5 *s_sc(3)^2/r^2 -3));

a_J2 = [a_j2_x; a_j2_y; a_j2_z];                         % [km/s^2]

% J3 Effect
kz = 1/2 * J3 * mi * R_pl^3 / r^5;

a_j3_x = kz*5*s_sc(1)/r*(7*(s_sc(3)/r)^3-3*(s_sc(3)/r));
a_j3_y = kz*5*s_sc(2)/r*(7*(s_sc(3)/r)^3-3*(s_sc(3)/r));
a_j3_z = kz*3*(1-10*(s_sc(3)/r)^2+35/3*(s_sc(3)/r)^4);

a_J3 = [a_j3_x; a_j3_y; a_j3_z];                        % [km/s^2]

% J4 Effect
kj = 5/8 * J4 * mi * R_pl^4 / r^6;
a_j4_x = kj*(3 - 42*(s_sc(3)/r)^2 + 63*(s_sc(3)/r)^4)*s_sc(1)/r;
a_j4_y = kj*(3 - 42*(s_sc(3)/r)^2 + 63*(s_sc(3)/r)^4)*s_sc(2)/r;
a_j4_z = kj*(15 - 70*(s_sc(3)/r)^2 + 63*(s_sc(3)/r)^4)*s_sc(1)/r;

a_J4 = [a_j4_x; a_j4_y; a_j4_z];                        % [km/s^2]

% rotation in ECI
a_J2 = B'*a_J2;
a_J3 = B'*a_J3;
a_J4 = B'*a_J4;

%% Out of gravity perturb
aJ2_car = a_J2 + a_J3 + a_J4;

%% ET capture thrust
mjd2000 = data.InitDay + t/86400; 

% uplanet
DEG2RAD=pi/180;


KM  = 149597870.66;

T   = (mjd2000 + 36525)/36525.00;
TT  = T*T;
TTT = T*TT;
KepS(1)=T*0;

KepS(1) = 1.5236883990;
KepS(2) = 0.093312900 + 0.0000920640*T - 0.0000000770*TT;
KepS(3) = 1.850333333333333330 - 6.75e-4*T + 1.26111111111111111e-5*TT;
KepS(4) = 4.87864416666666667e+1 + 7.70991666666666667e-1*T - 1.38888888888888889e-6*TT - 5.33333333333333333e-6*TTT;
KepS(5) = 2.85431761111111111e+2 + 1.069766666666666670*T +  1.3125e-4*TT + 4.13888888888888889e-6*TTT;
XM   = 1.91398585e+4 + 1.80805555555555556e-4*T + 1.19444444444444444e-6*TT;
KepS(6) = 3.19529425e2 + XM*T;

KepS(1)   = KepS(1)*KM;       % a [km]
KepS(3:6) = KepS(3:6)*DEG2RAD;    % Transform from deg to rad
KepS(6)   = mod(KepS(6),2*pi);
% XMU  = (XM*DEG2RAD/(864000*365250))^2*kep(1)^3;
phi    = KepS(6);          % phi is the eccentric anomaly, uses kep(6)=M as a first guess
     
for i=1:5
    g       = KepS(6)-(phi-KepS(2)*sin(phi)); 
	g_primo = (-1+KepS(2)*cos(phi));
 	phi     = phi-g/g_primo;   % Computes the eccentric anomaly kep
end 

theta=2*atan(sqrt((1+KepS(2))/(1-KepS(2)))*tan(phi/2));

KepS(6)=theta;

rM2S = - kep2car_r_only(KepS);         % retriving the sun position vector wrt the planet
light = los(R, rM2S, data.R_pl);       % checking the line of sight

if light
    [aT_car, thrust] = ETConstThrust(rM2S, R, V, M, data);            % thrust [N]
    th_range = ThRange(R, V, data);

    if not(th_range)
        aT_car = [0; 0; 0];
        thrust = 0;
    end

else
    aT_car = [0; 0; 0];
    thrust = 0;
end

%Isp = data.Isp_fun(thrust);

% Derivative of the states
dY(1:3) = V;
dY(4:6) = - data.mi/norm(R)^3.*R + aT_car + aJ2_car;
%dY(7) = - norm(thrust)/Isp/data.g0;
dY(7) = - norm(thrust)/4300/data.g0;


dY = dY';

end