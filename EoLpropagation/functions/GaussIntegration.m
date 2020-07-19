function [dY, R] = GaussIntegration(t, Y, data)
%{
CarIntegration.m - OdeFcn that computes the perturbed position and velocity
            vector using the direct cartesian approach. 
           
%}

% States
a = Y(1);
e = Y(2);
i = Y(3);
% RAAN = Y(4);
PA = Y(5);
TA = Y(6);

mi = data.mi;

% DayActual = data.InitDay + t/(24*3600);        % updated time [days]

[R,V] = kep2car (Y, mi);
r = norm(R);
v = norm(V);

%% Solar Radiation Pressure
if data.switchers.SRP
    a_SRP = SRP(t, R, data);
else
    a_SRP = [0; 0; 0];
end

%% Perturbation due to gravity
if data.switchers.grav
    [a_J2, a_J3, a_J4] = GravPerturb(t, R, data);
else
    a_J2 = [0; 0; 0];
    a_J3 = [0; 0; 0];
    a_J4 = [0; 0; 0];
end

%% Moon Perturbation
% if data.switchers.Moon
%     a_Phobos = MoonPerturb(DayActual, R, data, 'Phobos');
% %     a_Deimos = MoonPerturb(DayActual, R, mi_moon, 'Deimos');
% else
%     a_Phobos = [0; 0; 0];
% %     a_Deimos = [0; 0; 0];
% end


%% Pertubations in TNH frame
% Basis of TNH frame
tvers = V./v;
H = cross(R, V);
hvers = H./norm(H);
nvers = cross(hvers, tvers);

% Rotation to TNH frame
a_SRP_tnh = [tvers nvers hvers]'*a_SRP;
a_J2_tnh = [tvers nvers hvers]'*a_J2;
a_J3_tnh = [tvers nvers hvers]'*a_J3;
a_J4_tnh = [tvers nvers hvers]'*a_J4;
% a_Phob_tnh = [tvers nvers hvers]'*a_Phobos;

a_tnh = a_SRP_tnh + a_J2_tnh + a_J3_tnh + a_J4_tnh; % +  a_Phob_tnh;

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

%% State derivatives
dY(1) = da;
dY(2) = de;
dY(3) = di;
dY(4) = dRAAN;
dY(5) = dPA;
dY(6) = dTA;

%% Extra outputs
% parout{1} = R;

dY = dY';


end