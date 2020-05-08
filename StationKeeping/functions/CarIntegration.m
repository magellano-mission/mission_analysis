function [dY, parout] = CarIntegration(t, Y, data)
%{
CarIntegration.m - OdeFcn that computes the perturbed position and velocity
            vector using the direct cartesian approach. 
           
%}

% States
R = Y(1:3);
V = Y(4:6);

% input checks
if isrow(R)
    R = R';
end
if isrow(V)
    V = V';
end

mi = data.mi;

DayActual = data.InitDay + t/(24*3600);        % updated time [days]

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
if data.switchers.Moon
    a_Phobos = MoonPerturb(DayActual, R, data, 'Phobos');
%     a_Deimos = MoonPerturb(DayActual, R, mi_moon, 'Deimos');
else
    a_Phobos = [0; 0; 0];
%     a_Deimos = [0; 0; 0];
end

%% total perturbation
a_p = a_SRP + a_J2 + a_J3 + a_J4 + a_Phobos; % + a_Deimos;

% State derivatives
dY(1:3) = V;
dY(4:6) = - mi/norm(R)^3.*R + a_p;

dY = dY';

%% parout
[kep] = car2kep(R, V, mi);
parout{1} = kep(1);
parout{2} = kep(2);
parout{3} = kep(3)*180/pi;
parout{4} = kep(4)*180/pi;
parout{5} = kep(5)*180/pi;
parout{6} = kep(6)*180/pi;


end
