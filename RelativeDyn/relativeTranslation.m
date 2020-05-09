function[dX]=relativeTranslation(t,X,data)
%(!!!!! Everything is in m )
% Extract Parameters 

a  = data.a ; % semi-major axis of leader orbit % [m]
e  = data.e ; % eccentricity of leader orbit 
mu = data. mu; % gravitational parameter 

dX = zeros(length(X),1);  % Initialize state vector 
% Extract State Parameters 
x      = X(1);
y      = X(2);
z      = X(3);
xdot   = X(4);
ydot   = X(5);

thL    = X(7);
                                                                                                            
thLdot = sqrt(mu/(a^3*(1-e^2)^3)) * (1+e*cos(thL))^2;

rL    = a*(1-e^2)/(1+e*cos(thL));        % Leader Position in inertial frame 

h     = rL^2*thLdot;
rLdot = mu/h*e*sin(thL);

thLddot = - 2*rLdot*thLdot/rL;


dX(1:3) = X(4:6);

dX(4) =  2*thLdot*ydot + thLddot*y + thLdot^2*x - mu*(rL+x)/((rL+x)^2+y^2+z^2)^(3/2) + mu/rL^2; 

dX(5) = -2*thLdot*xdot - thLddot*x + thLdot^2*y - mu* y/((rL+x)^2+y^2+z^2)^(3/2);

dX(6) = - mu * z / ((rL+x)^2+y^2+z^2)^(3/2);

dX(7) = thLdot; %thetadot 








