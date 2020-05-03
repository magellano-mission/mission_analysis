function [SAA_] = SAA_fun(a,e,i,om,OM,th,In_time,mu,n)

% Orbital period for circular orbit
T = 2*pi*sqrt( a^3 / mu );

[r_0, v_0] = kep2car2(a, e, i, om, OM, th, mu);  

y01 = [r_0; v_0];

[~, y1] = ode113( @(t, y) elliptical_orbit(t, y, mu), linspace(In_time,T+In_time,500), y01, odeset('RelTol', 1e-13, 'AbsTol', 1e-14));

t_SAA = linspace(In_time,T+In_time,500);

SAA_ = zeros(length(t_SAA),1);

for i = 1:length(t_SAA)
    
    t = t_SAA(i);
    
    [k, mu_s] = uplanet(t/3600/24, n);
    
    [r, ~] = kep2car2(k(1), k(2), k(3), k(4), k(5), k(6), mu_s); 
    
    SAA_(i,1) =  acos(dot(y1(i,1:3),(r'+y1(i,1:3))) / (norm((r'+y1(i,1:3)))*norm(y1(i,1:3)))) ;
    
end

figure

plot(t_SAA,rad2deg(SAA_),'Linewidth', 2, 'DisplayName', 'Transfer arc', 'Color', 'r')
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 11)
xlabel('Time [day]')
ylabel('Angle [deg]')
title(' Planet-s/c-Sun angle over one period')
grid on; grid minor;

end

