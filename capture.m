function [dt_flyby, dv_tot] = capture(kep_cap, delta, parameters, mu)
% function [rp,dt_flyby, dv_tot] = capture(ID1,ID2,ID3,R_atm,flyby_date,dep_date,arr_date, parameters, mu)
%PLOT_GA 
%Plot of asymptotes and hyperbolic path in occurrence of the flyby
%
% INPUT
%   - IDx [-] Planets' codes
%   - R_atm [km] Atmospheric upper bound (Constraint)
%   - x_date [days] MJD2000 date of flyby/departure/arrival 
%OUTPUT
%   - rp [km] minimum radius of approach of hyperbolic path (from center)
%   - dt_flyby [hours] TOF across Flyby planet's SOI
%   - dv_tot [km/s] Total deltaV provided by powered g.a. (dv@peric.+ g.a.)
%
%
% Authors : Apeksha Veeranna Roopashree / Sergio Bonnacorsi /
%                     Alberto Chiaradia / Louis Aucouturier
%

%Discretization of hyperbola
n_iter = 1000;

%Computation of the two hyperbolas

e_m = 1/sin(delta/2);

h_m = sqrt(mu*rp*(1+e_m));

a_m =((h_m^2)/mu)/(e_m^2-1);

theta_infm = acos(1/e_m);

%SOI definition
rSOI = nrM * (mu/muS)^(2/5);

correc_m_f = @(alpha) correc_SOI_f(rSOI,alpha,a_m,e_m,theta_infm,mu);

correc_m = fminsearch(correc_m_f,0);

ltheta_m = linspace(0,theta_infm+correc_m,n_iter);

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
dt_flyby = (dt_m)/3600;

%Hyperbolas plot
l_pos_Vect_m = zeros(3,n_iter);

kep_hyp = kep_cap;
kep_hyp(1) = a_m;
kep_hyp(2) = e_m;

for i = 1:n_iter
    kep_hyp(6) = ltheta_m(i);
    [r_m,~]= kep2car(kep_hyp,mu);
    l_pos_Vect_m(:,i) = r_m;
end

figure
axis equal
hold all

r = rSOI;
xc = 0;
yc = 0;

theta = linspace(0,2*pi);
x = r*cos(theta) + xc;
y = r*sin(theta) + yc;

plot(x,y)
plot3(l_pos_Vect_m(1,:),l_pos_Vect_m(2,:), l_pos_Vect_m(3,:)), hold on

hold off
end

% Function used to get a correction on the theta for the hyperbola arc to
% reach the Sphere of Influence
function val = correc_SOI_f(rSOI,alpha,a,e,theta,muM)
    [r,~]= kep2car(a,e,0,0,0,theta+alpha,muM);
    val = (norm(r)-rSOI)^2;
end



