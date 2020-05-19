function[rr,vv]=kep2car2(kep,mu)
% Transformation from Keplerian elements to cartesian coordinates
%
% -------------------------------------------------------------------------
% Input arguments:
% a     [1x1] semi-major axis           [km]
% e     [1x1] eccentricity              [-]
% i     [1x1] inclination               [rad]
% OM    [1x1] RAAN                      [rad]
% om    [1x1] argument of periapsis     [rad]
% th    [1x1] true anomaly              [rad]
% mu    [1x1] gravitational parameter   [km^3/s^2]
%
% -------------------------------------------------------------------------
% Output arguments:
% rr [3x1] position vector [km]
% vv [3x1] velocity vector [km/s]

a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
th = kep(6);


R3_OM=[cos(OM) sin(OM) 0
      -sin(OM) cos(OM) 0
          0       0    1];
R1_i= [1    0       0
       0   cos(i) sin(i)
       0  -sin(i) cos(i)];
R3_om=[cos(om) sin(om) 0
      -sin(om) cos(om) 0
          0       0    1];
R=[R3_om*R1_i*R3_OM]';
if e<=1
    p=a*(1-e^2);
else
    p=a*(e^2-1);
end
r=p/(1+e*cos(th));

r_pf=r*[cos(th)
        sin(th)
           0   ];
v_pf=sqrt(mu/p)*[ -sin(th)
                 e+cos(th)
                      0   ];

rr=R*r_pf;
vv=R*v_pf;

end