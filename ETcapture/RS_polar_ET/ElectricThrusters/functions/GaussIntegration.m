function [dY, R] = GaussIntegration(t, Y, data)
%{
CarIntegration.m - OdeFcn that computes the perturbed position and velocity
            vector using the direct cartesian approach.           
%}
t
% States
a = Y(1);
e = Y(2);
i = Y(3);
% RAAN = Y(4);
PA = Y(5);
TA = Y(6);
M = Y(7);

mi = data.mi;

[R,V] = kep2car(Y, mi);
r = norm(R);
v = norm(V);

T = data.T;
if data.direction == "tangential"
    aT_tnh = [T; 0; 0]/M*1e-3;             
elseif data.direction == "perpendicular"
    aT_tnh = [0; 0; data.T]/M*1e-3;
end

%% Perturbation due to gravity
[a_J2, ~, ~] = GravPerturb(t, R, data);

%% Pertubations in TNH frame
% Basis of TNH frame
tvers = V./v;
H = cross(R, V);
hvers = H./norm(H);
nvers = cross(hvers, tvers);

% Rotation to TNH frame
a_J2_tnh = [tvers nvers hvers]'*a_J2;

a_tnh = a_J2_tnh + aT_tnh;

% Components of accelerations in TNH frame
at = a_tnh(1);
an = a_tnh(2);
ah = a_tnh(3);

% Parameters for the Gauss Equations
n = sqrt(mi/(a^3));
b = a*sqrt(1 - e^2);
h = n*a*b;
TA_s = TA + PA;

%% Gauss equations
da = 2*a^2*v/mi * at;
de = 1/v*(2*(e+cos(TA))* at - r/a*sin(TA) * an);
di = r* cos(TA_s)/h * ah;
dRAAN = r*sin(TA_s)/(h*sin(i)) * ah;
dPA = 1/(e*v)*(2*sin(TA) * at + (2*e+r/a*cos(TA)) * an) - r*sin(TA_s)*cos(i)/(h*sin(i)) * ah;
dTA = h/(r^2) - 1/(e*v)*(2*sin(TA)*at + (2*e + r/a*cos(TA))*an);
dM = - norm(data.T)/data.Isp/data.g0;

%% State derivatives
dY(1) = da;
dY(2) = de;
dY(3) = di;
dY(4) = dRAAN;
dY(5) = dPA;
dY(6) = dTA;
dY(7) = dM;

% Output
dY = dY';



end


