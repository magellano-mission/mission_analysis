function[dX,parout]=coupledRD(t,X,data)
%(!!!!! Everything is in m )
% Extract Parameters 

a  = data.a ; % semi-major axis of leader orbit % [km]
e  = data.e ; % eccentricity of leader orbit 
i  = data.i; 
OM = data.OM; 
om = data.om; 
th = X(14); 
mu = data. mu; % gravitational parameter 
JT = diag(data.JT); JTinv = diag(1./data.JT);
JL = diag(data.JT); JLinv = diag(1./data.JL);
dX = zeros(length(X),1);  % Initialize state vector 




% Extract State Parameters 
% Relative Position 
x      = X(1);
y      = X(2);
z      = X(3);
% Relative velocity 
xdot   = X(4);
ydot   = X(5);

% Relative angular velocity 
w1     = X(7); 
w2     = X(8);
w3     = X(9);

% Relative quaternion (q0,q1,q2,q3)
quat = [X(10);X(11);X(12);X(13)];
q0     = quat(1);
q1     = quat(2);
q2     = quat(3);
q3     = quat(4);

% Leader true anomaly 
thL    = X(14);

%%%%%%%%%%%%%%% Preliminary Computations

% Position Vector ( Cartesian Coordinates)
[rr,vv]=kep2car([a e i OM om th],mu);
hv = cross(rr,vv); 
rL = norm(rr); 
h  = norm(hv); 
% Leader Position, True Anomaly , leader frame angular velocity 
thLdot = sqrt(mu/(a^3*(1-e^2)^3)) * (1+e*cos(thL))^2;
rLdot = mu/h*e*sin(thL);
thLddot = - 2*rLdot*thLdot/rL;

% Attitude Parametrization 
A = [ q1^2-q2^2-q3^2+q0^2              2*(q1*q2-q3*q0)            2*(q1*q3+q2*q0)
       2*(q1*q2+q3*q0)             -q1^2+q2^2-q3^2+q0^2           2*(q2*q3-q1*q0) 
       2*(q1*q3-q2*q0)                2*(q2*q3+q1*q0)         -q1^2-q2^2+q3^2+q0^2];
 

om  =  [w1;w2;w3];
omL =  [0;0;thLdot];
Q = [ -q1 -q2 -q3
       q0 -q3  q2
       q3  q0 -q1 
      -q2  q1  q0 ]; 
% Relative Position and Velocity Vector 

dX(1:3) = X(4:6);

dX(4) =  2*thLdot*ydot + thLddot*y + thLdot^2*x - mu*(rL+x)/((rL+x)^2+y^2+z^2)^(3/2) + mu/rL^2; 

dX(5) = -2*thLdot*xdot - thLddot*x + thLdot^2*y - mu* y/((rL+x)^2+y^2+z^2)^(3/2);

dX(6) = - mu * z / ((rL+x)^2+y^2+z^2)^(3/2);

% Orbital True Anomaly 
dX(14) = thLdot;  

% Rotational  Dynamics 
dX(7:9) = JLinv*( JL*A*JTinv*(cross( -A'*(om+omL) , JT*(A'*(om+omL)))) -cross(JL*omL,om) +(cross(omL,JL*omL))) ;

% Quaternion Kinematics
dX(10:13) = 1/2*Q*A'*om; 



% Feature Points dynamics 
% rho_i = rho_cm + cross(om,Pi) 
% Pi : posistion vector from target Cm to feature point, expressed in
% leader frame 

nF = length(X(15:end))/3; 
parout = zeros(1,nF*3);
for i = 1 :nF  
    dX(14+3*i-2:14+3*i) = cross(om,X(14+3*i-2:14+3*i));
    parout(3*i-2:3*i) = X(1:3)+ X(14+3*i-2:14+3*i);
end 





