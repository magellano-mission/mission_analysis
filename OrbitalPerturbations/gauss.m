function[dx, parout]=gauss(t,x)

% INPUT
% t   [1]   : epoch  of the integration
% x [1x6] : State vector of the keplerian parameters [a e i OM om theta]
%                    - a = kep(1) [km]
%                    - e = kep(2) [-]
%                    - i = kep(3) [rad]
%                    - OM = kep(4) [rad]
%                    - om = kep(5) [rad]
%                    - theta = kep(6) [rad]
% x_id [1x6] : State vector of the ideal keplerian parameters [a e i OM om theta]
%
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
J3 = 3.1450e-5; %Mars J3
J4 = -1.5377e-5; %Mars J4

%Preliminary calculations
theta_star = theta+om;
n = sqrt(mu_m/a^3);
b = a*sqrt(1-e^2);
h = n*a*b;
p = b^2/a;
r = p/(1+e*cos(theta));
v = sqrt( 2*mu_m/r  - mu_m/a);

% Cartesian Coordinates of the satellite
[rr, vv] = kep2car(kep, mu_m);
[altitude, R_real] = MCI2altitude(rr); %[km]

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

% J3 Effect
kz = 1/2 * J3 * mu_m * R_real^3 / r^5;

a_j3_x = kz*5*rr(1)/r*(7*(rr(3)/r)^3-3*(rr(3)/r));
a_j3_y = kz*5*rr(2)/r*(7*(rr(3)/r)^3-3*(rr(3)/r));
a_j3_z = kz*3*(1-10*(rr(3)/r)^2+35/3*(rr(3)/r)^4);

a_J3 = [a_j3_x; a_j3_y; a_j3_z]; % [km/s^2]

% J4 Effect
kj = 5/8 * J4 * mu_m * R_real^4 / r^6;
a_j4_x = kj*(3 - 42*(rr(3)/r)^2 + 63*(rr(3)/r)^4)*rr(1)/r;
a_j4_y = kj*(3 - 42*(rr(3)/r)^2 + 63*(rr(3)/r)^4)*rr(2)/r;
a_j4_z = kj*(15 - 70*(rr(3)/r)^2 + 63*(rr(3)/r)^4)*rr(1)/r;

a_J4 = [a_j4_x; a_j4_y; a_j4_z]; % [km/s^2]

% Phobos data
M_phobos = 1.07e16; %Phobos mass
G = astroConstants(1);
mu_phobos = M_phobos * G;

% Retrieving moon position
%mjd = date2mjd2000(data.date) + t/day2sec;
stat_phobos = [9376, 0.0151, deg2rad(1.075) deg2rad(317.671) deg2rad(150.057) 0];
t_span = [0 t+1]; %the first one should be the time (in JD) in which we know position
[~, x_phobos] = ode113(@gauss_sat, t_span,  stat_phobos);
[rr_phobos, ~]=kep2car(x_phobos(end,:),mu_m);
r_phobos = norm(rr_phobos); %[km]

% Retrieving relative position moon / spacecraft
rr_rel_m = rr_phobos - rr;
r_rel_m = norm(rr_rel_m);

% Phobos gravitational acceleration
a_phobos = mu_phobos * (rr_rel_m/r_rel_m^3 - rr_phobos/r_phobos^3);

% a_phobos = 0;

% Deimos data
M_deimos = 1.4762e15;
mu_deimos = M_deimos * G;

% Retrieving moon position
%mjd = date2mjd2000(data.date) + t/day2sec;
stat_deimos = [23458, 0.0002, deg2rad(1.788) deg2rad(316.657) deg2rad(260.729) 0];
t_span = [0 t+1]; %the first one should be the time (in JD) in which we know position
[~, x_deimos] = ode113(@gauss_sat, t_span,  stat_deimos);
[rr_deimos, ~]=kep2car(x_deimos(end,:),mu_m);
r_deimos = norm(rr_deimos); %[km]

% Retrieving relative position moon / spacecraft
rr_rel_d = rr_deimos - rr;
r_rel_d = norm(rr_rel_d);

% Deimos gravitational acceleration
a_deimos = mu_deimos * (rr_rel_d/r_rel_d^3 - rr_deimos/r_deimos^3);

% a_deimos = 0;

% SRP acceleration

ibody = 4; %Mars
[kep,ksun] = uplanet(t/(24*3600), ibody); % date in days
[rr_sun,~] = kep2car (kep,ksun);
r_sun = norm(rr_sun);

% identify s/c-sun vector

rr_rel_sun = -(rr_sun+rr);
R_sc_sun = norm(rr_rel_sun);

% SRP perturbation

R_AU = 1.52*1.496*10^8;
c = 3*10^8;
W =  609.5000; %average value between min and max
p_sr = W/c;
c_r = 1; %reflection coefficient
A_m = 5; % Area to mass ratio

%Asrp = 0.001*p_sr*(R_AU^2)*c_r*A_m/R_sc_sun^2;
Asrp = p_sr*(R_AU^2)*c_r*A_m/R_sc_sun^2*10^-6; %check order of magnitude
a_srpi = -Asrp*rr_rel_sun(1)/R_sc_sun;
a_srpj = -Asrp*rr_rel_sun(2)/R_sc_sun;
a_srpk = -Asrp*rr_rel_sun(3)/R_sc_sun;

a_srp = [a_srpi; a_srpj; a_srpk];


% Sun perturbation

rr_rel_sun = rr_sun - rr;
r_rel_sun = norm(rr_rel_sun);

% Deimos gravitational acceleration
a_sun = ksun * (rr_rel_sun/r_rel_sun^3 - rr_sun/r_sun^3);

% DRAG perturbation

wMars = 2*pi/(24.6*60*60);
v_rel = vv - cross([0;0;wMars],(rr-R_real));
V_REL = norm(v_rel);
Ac_m = 1; % cross area to mass ratio
cD = 1; % estimated drag coefficient

if altitude <50
    rho = 0.01442*exp(-0.08819*altitude);
else
    rho = 0.001039*exp(-0.04834*altitude);
end
if altitude>1000
%     rho = 0.001039*exp(-0.04834*1000);
    rho = 0;
end

a_drag = 0.5*rho*V_REL*v_rel*Ac_m*cD;

% TOTAL ACCELERATIONS
a_tot = a_J2 + a_J3 + a_J4 + a_phobos + a_deimos + a_srp + a_drag + a_sun;
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

a_primary = mu_m/norm(rr)^2;
parout = [a_primary, norm(a_J2), norm(a_J3), norm(a_J4), norm(a_phobos), norm(a_deimos), norm(a_sun), norm(a_srp), norm(a_drag)];
end
