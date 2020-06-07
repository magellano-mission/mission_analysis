function[dx]=gauss_sat(t,x)

% gauss for the two satellites. No perturbation is accounted

% INPUT
% t   [1]   : epoch  of the integration
% kep [1x6] : State vector of the keplerian parameters [a e i OM om theta]
%                    - a = kep(1) [km]
%                    - e = kep(2) ,
%                    - i = kep(3) [rad]
%                    - om = kep(4) [rad]
%                    - OM = kep(5) [rad]
%                    - theta = kep(6) [rad]
% OUTPUT :
% dkepdt[6x1]: Derivative of the state (transposed to match the dimensions
%              required by the numerical integrator)
%
%
% Other m-files required: astroConstants
%                         ECI2altitude
%                         kep2car 

[n,m]=size(x);
dx=zeros(n,m);

% Keplerian Parameter Propagation

% Extract the keplerian parameters from the state vector
a=x(1);
e=x(2);
i=x(3);
om=x(5);
theta=x(6);
kep = x(1:6);

%Physical data and Time definition
mu_m = astroConstants(14);
J2 = 0.00196; %Mars J2

%Preliminary calculations
theta_star = theta+om;
n = sqrt(mu_m/a^3);
b = a*sqrt(1-e^2);
h = n*a*b;
p = b^2/a;
r = p/(1+e*cos(theta));
v = sqrt( 2*mu_m/r  - mu_m/a);

% Cartesian Coordinates of the satellite
[rr, vv] = kep2car2(kep, mu_m);
[~, R_real] = MCI2altitude(rr); %[km] needed for?

% Rotation Matrix
t_hat = vv/norm(vv);
h_hat = cross(rr,vv)/norm(cross(rr,vv));
n_hat = cross(h_hat,t_hat);
A = [ t_hat  n_hat  h_hat ];

% J2 Effect
k = 3/2 * J2 * mu_m * R_real^2 / r^4;

a_j2_x = k*(rr(1)/r *(5 *rr(3)^2/r^2 -1));
a_j2_y = k*(rr(2)/r *(5 *rr(3)^2/r^2 -1));
a_j2_z = k*(rr(3)/r *(5 *rr(3)^2/r^2 -3));

a_J2 = [a_j2_x; a_j2_y; a_j2_z]; % [km/s^2]

% Conversion day to seconds
%day2sec = 86400;

a_tot = a_J2;
a_tnh = A' * a_tot;

at = a_tnh(1);  % Thrust aligned with velocity vector (converted to km/s^2)
an = a_tnh(2);
ah = a_tnh(3);

dx(1) = 2*a^2*v*at/mu_m;
dx(2) = 1/v*(2*(e+cos(theta))*at-r/a*sin(theta)*an);
dx(3) = (r*cos(theta_star)/h)*ah;
dx(4) = r*sin(theta_star)*ah/(h*sin(i));
dx(5) = 1/(e*v)*(2*sin(theta)*at + (2*e+r*cos(theta)/a)*an) - r*sin(theta_star)*cos(i)*ah/(h*sin(i));
dx(6) = h/r^2 - 1/(e*v)*(2*sin(theta)*at + (2*e+r*cos(theta)/a)*an );
end
