function dv = SKmaneuver(kep1, data)

mi = data.mi;
a1 = kep1(1);
e1 = kep1(2);
i1 = kep1(3)*pi/180;
% a1 = 10600;
% e1 = 0.1;
% i1 = 26*pi/180;
th1 = pi;               % the maneuver should be done at pericenter for low dv

rcirc = data.sma;
icirc = data.inc;
% rcirc = 10500;
% icirc = 25*pi/180;

%%%% first dv: change of plane and transfer orbit
p1 = a1*(1 - e1^2);
ra1 = p1/(1-e1);
v1 = sqrt(2*mi*(1/ra1 - 1/(2*a1)));

ht = sqrt(2*mi*ra1*rcirc/(ra1+rcirc));
v1t = ht/ra1;

dv1 = sqrt(v1^2 + v1t^2 - 2*v1*v1t*cos(abs(i1 - icirc)));

%%%% first dv: from transfer orbit to circular
v2t = ht/rcirc;

v2 = sqrt(mi/rcirc);

dv2 = abs(v2 - v2t);

%%%% total dv
dv = dv1 + dv2;

