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

date = data.initialDate;
mi_earth = data.mi;

DayInit = date2mjd2000(date);             % initial time [days]
DayActual = DayInit + t/(24*3600);        % updated time [days]

%% Solar Radiation Pressure
if data.switchers.SRP
    a_SRP = SRP(t, R, data);
else
    a_SRP = [0; 0; 0];
end

%% Perturbation due to J2
if data.switchers.J2
    [a_J2, a_J3, a_J4] = GravPerturb(t, R, data);
else
    a_J2 = [0; 0; 0];
end

%% Moon Perturbation
if data.switchers.Moon
    a_Phobos = MoonPerturb(DayActual, R, mi_moon, 'Phobos');
    a_Deimos = MoonPerturb(DayActual, R, mi_moon, 'Deimos');
else
    a_Phobos = [0; 0; 0];
    a_Deimos = [0; 0; 0];
end

%% total perturbation
a_p = a_SRP + a_J2 + a_J3 + a_J4 + a_Phobos + a_Deimos;

%% State derivatives
dY(1:3) = V;
dY(4:6) = - mi_earth/norm(R)^3.*R + a_p;

dY = dY';

%% parout
[orbital_parameters] = car2kep(R, V, mi_earth);
parout{1} = orbital_parameters.a;
parout{2} = orbital_parameters.e;
parout{3} = orbital_parameters.i;
parout{4} = orbital_parameters.RAAN;
parout{5} = orbital_parameters.PA;
parout{6} = orbital_parameters.TA;


end
