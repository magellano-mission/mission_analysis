function [EAA_] = Earth_visibility_angle(Earth_time, Mars_time)

[k_E, mu_s] = uplanet(Earth_time, 3);
[r_E, ~] = kep2car(k_E(1), k_E(2), k_E(3), k_E(4), k_E(5), k_E(6), mu_s);
        
[k_E, ~] = uplanet(Mars_time, 4);
[r_M, ~] = kep2car(k_E(1), k_E(2), k_E(3), k_E(4), k_E(5), k_E(6), mu_s);
           
[~, ~, ~, ~, VI, ~, ~, ~] = lambertMR(r_E, r_M, (Mars_time - Earth_time)*24*3600 , mu_s);

y01 = [r_E; VI'];

[~, y1] = ode113( @(t, y) elliptical_orbit(t, y, mu_s) , linspace(Earth_time, Mars_time, 500)*24*3600, y01, odeset('RelTol', 1e-13, 'AbsTol', 1e-14));

t_E_EAA = linspace(Earth_time, Mars_time, 500);

for i = 1:length(t_E_EAA)
    
    tE = t_E_EAA(i);
    
    [k_E, ~] = uplanet(tE, 3);    
    [r, ~] = kep2car(k_E(1), k_E(2), k_E(3), k_E(4), k_E(5), k_E(6), mu_s);    
    rr_E_EAA(:, i) = r;
    
end

vec = y1(:,1:3) - rr_E_EAA';

EAA_ = zeros(length(t_E_EAA),1);

for i=1:length(t_E_EAA)
    
EAA_(i,1) =  acos(dot(y1(i,4:6),-vec(i,:)) / (norm(-vec(i,:))*norm(y1(i,4:6)))) ;

end

figure
plot(t_E_EAA,rad2deg(EAA_),'Linewidth', 2, 'DisplayName', 'Transfer arc', 'Color', 'r')
grid on; grid minor;
xlabel('Time[day]');
ylabel('Angle [deg]');
title('Earth-s/c-\underline v ')

end