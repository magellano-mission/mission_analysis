function rr = refplan2car( hh, th )
% Transformation from reference plane to cartesian

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

r_pf=r*[cos(th)
        sin(th)
           0   ];
rr=R*r_pf;

end