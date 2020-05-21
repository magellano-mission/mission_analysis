function [rr, vv] = refplane2car( r, z, r1vers, vr, vt, hh, th, muS )
% Transformation from reference plane to cartesian
ee = r1vers;

r_pf=[r*cos(th)
      r*sin(th)
           z   ];
v_pf = [vr
        vt
        0];
%Calcolo dell'inclinazione.
i = acos(hh(3)/norm(hh)); %versor

%Calcolo della linea dei nodi.
nn = cross([0,0,1]',hh);
n = norm(nn);

%Calcolo dell'ascensione retta del nodo ascendente.
OM = acos(nn(1)/n);
if nn(2)<0
    OM = 2*pi - OM;
end

%Calcolo dell'anomalia del pericentro.
om = acos(dot(nn/n,ee));
if ee(3)<0
    om = 2*pi-om;
end
om = 0;

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


rr=R*r_pf;
vv=R*v_pf;
end