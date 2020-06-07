function [rr, vv]=kep2car(stat, mi)
% [rr,vv] = kep2car(a,e,i,Om,om,th,mi) 
% Computation of the cartesian position and velocity vectors starting from
% keplerian parameters.
% 
% INPUTS
% stat  : [1x6]  vector of keplerian parameters of the orbit [a e i Om om theta]
% mu    : [1]    gravitational parameter of the main body   
% OUTPUTS
% rr    : [1x3]  position vector
% vv    : [1x3]  velocity vector
%
% 

%Extracting parameters from stat
a=stat(1); e=stat(2);
i=stat(3); Om=stat(4); 
om=stat(5); th=stat(6);

%Perifocal frame computation PF
p = a*(1-e^2);
if e > 1
    p = a*(e^2 - 1);
end

%Magnitude of the position vector
r=p/(1+e*cos(th));

%rr and vv in PF
rr_pf= r*[cos(th); sin(th); 0];
vv_pf=(sqrt(mi/p))*[-sin(th); e+cos(th); 0];

% ROTATION MARIX
%Rotation matrixes for Om, om and i
R_Om=[cos(Om) sin(Om) 0; -sin(Om) cos(Om) 0; 0 0 1];
R_i=[1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R_om=[cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

%From PF to GE

R=R_Om'*R_i'*R_om';

% Computation of rr and vv in GE
rr=R*rr_pf;
vv=R*vv_pf;


end