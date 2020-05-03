function [dt_hyp, dv_req, theta_inf, kep_hyp] = capture(kep_cap, delta,rM, mu)
%Plot of asymptotes and hyperbolic path in occurrence of the flyby
%   INPUT
%   - delta: deflection angle (free parameters JUST TO RUN SIMULATIONS NOW)
%   - rM [km] planet position wrt SUN
%   - kep_cap keplerian parameters of the desired parking orbit
%OUTPUT
%   - rp [km] minimum radius of approach of hyperbolic path (from center)
%   - dt_flyby [hours] TOF across Flyby planet's SOI
%   - dv_req [km/s] impulsive delta_v

muS = astroConstants(4);
%Discretization of hyperbola
n_iter = 1000;

%Computation of the hyperbola

rp = kep_cap(1)*(1 - kep_cap(2));

e_m = 1/sin(delta/2);

h_m = sqrt(mu*rp*(1+e_m));

a_m =((h_m^2)/mu)/(e_m^2-1);

theta_inf = acos(1/e_m);

nrM = norm(rM); 
%SOI definition
rSOI = nrM * (mu/muS)^(2/5);

correc_m_f = @(alpha) correc_SOI_f(rSOI,alpha,a_m,e_m,theta_inf,mu);
correc_m = fminsearch(correc_m_f,0);
ltheta_m = linspace(0,theta_inf+correc_m,n_iter);

lt_m = zeros(1,n_iter);

%TOF computation
for i = 1:n_iter
    %eccentric anomalies
    f_m = 2*atanh(sqrt((e_m-1)/(e_m+1))*tan(ltheta_m(i)/2));
    m_m = e_m*sinh(f_m) - f_m;
    lt_m(i) = m_m*(h_m^3/mu^2)/((e_m^2-1)^(3/2));
end

%Computation of TOF
dt_m = abs(lt_m(1) - lt_m(n_iter));
dt_hyp = (dt_m)/3600;

%Hyperbolas plot
l_pos_Vect_m = zeros(3,n_iter);

kep_hyp = kep_cap;
kep_hyp(1) = a_m;
kep_hyp(2) = e_m;

[~, v_ell] = kep2car2(kep_cap, mu);
[~, v_hyp] = kep2car2(kep_hyp, mu);
dv_req = v_ell - v_hyp;

for i = 1:n_iter
    kep_hyp(6) = ltheta_m(i);
    [r_m,~]= kep2car2(kep_hyp,mu);
    l_pos_Vect_m(:,i) = r_m;
end

figure
hold all

r = rSOI;
xc = 0;
yc = 0;

theta = linspace(0,2*pi);
x = r*cos(theta) + xc;
y = r*sin(theta) + yc;

plot(x,y)
plot3(l_pos_Vect_m(1,:),l_pos_Vect_m(2,:), l_pos_Vect_m(3,:)), hold on


I = imread('Mars.jpg');                            % Mars image
RI = imref2d(size(I));
RI.XWorldLimits = [-180 180];                       % Mars image x sizes
RI.YWorldLimits = [-90 90];                         % Mars image y sizes
        rM = almanac('Mars','Radius','kilometers','sphere');
        [X, Y, Z] = ellipsoid(0, 0, 0, rM, rM, rM, 100); % spheric centered Mars
        planet = surf(X, Y, -Z,'Edgecolor', 'none');
        hold on
        set(planet,'FaceColor','texturemap','Cdata',I)
        axis equal


end

% Function used to get a correction on the theta for the hyperbola arc to
% reach the Sphere of Influence
function val = correc_SOI_f(rSOI,alpha,a,e,theta,muM)
    [r,~]= kep2car2([a,e,0,0,0,theta+alpha],muM);
    val = (norm(r)-rSOI)^2;
end



