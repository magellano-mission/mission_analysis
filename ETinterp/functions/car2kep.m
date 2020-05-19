function K = car2kep( rr, vv, mu )
% Transformation from cartesian coordinates to Keplerian elements
%
% [ a, e, i, OM, om, th] = car2kep( rr, vv, mu )
%
% -------------------------------------------------------------------------
% Input arguments:
% rr [3x1] position vector [km]
% vv[3x1] velocity vector [km/s]
% mu [1x1] gravitational parameter [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% a [1x1] semi-major axis [km]
% e [1x1] eccentricity [-]
% i[1x1] inclination [rad]
% OM [1x1] RAAN [rad]
% om [1x1] argument of periapsis [rad]
% th[1x1] true anomaly [rad]

%Calcolo il modulo di r e v.
r = norm(rr);
v = norm(vv);

%Calcolo dell'energia e del semiasse a.
E = (v^2)/2 -mu/r;
a = -mu/(2*E);

if a < 0
    a = abs(a);
end

%Calcolo dell'eccentricità.
ee = ((v^2-mu/r)*rr-(dot(rr,vv)*vv))/mu;
e = norm(ee);

%Calcolo del momento della quantità di moto.
hh = cross(rr,vv);
h = norm(hh);

%Calcolo dell'inclinazione.
i = acos(hh(3)/h);

%Calcolo della linea dei nodi.
nn = cross([0,0,1],hh);
n = norm(nn);

%Calcolo dell'ascensione retta del nodo ascendente.
OM = acos(nn(1)/n);
if nn(2)<0
    OM = 2*pi - OM;
end


%Calcolo dell'anomalia del pericentro.
om = acos(dot(nn,ee)/(n*e));
if ee(3) < 0
    om = 2*pi-om;
end

%Calcolo dell'anomalia vera.
vr = dot(vv,rr)/r;
th = acos( max( min( dot(ee, rr)/(e*r) , 1) , -1) );

if vr < 0 && e > 1
    th = - th;
else 
    th = 2*pi - th;
end

K = [a, e, i, OM, om, th];

end