function [EAA] = SAA_fun_tran(Earth_time, Mars_time)

[k_E, mu_s] = uplanet(Earth_time, 3);
[r_E, ~] = kep2car2([k_E(1), k_E(2), k_E(3), k_E(4), k_E(5), k_E(6)], mu_s);
        
[k_M, ~] = uplanet(Mars_time, 4);
[r_M, ~] = kep2car2([k_M(1), k_M(2), k_M(3), k_M(4), k_M(5), k_M(6)], mu_s);
           
[~, ~, ~, ~, VI, ~, ~, ~] = lambertMR(r_E, r_M, (Mars_time - Earth_time)*24*3600 , mu_s);

y01 = [r_E; VI'];

[~, y1] = ode113( @(t, y) elliptical_orbit(t, y, mu_s) , linspace(Earth_time, Mars_time, 500)*24*3600, y01, odeset('RelTol', 1e-13, 'AbsTol', 1e-14));

t_IT = linspace(Earth_time, Mars_time, 500);

rr_pl = zeros(3,length(t_IT));

EAA = zeros(length(t_IT),1);

for i = 1:length(t_IT)
    
    tE = t_IT(i);  
    [k_E, ~] = uplanet(tE, 3); 
    [r, ~] = kep2car2([k_E(1), k_E(2), k_E(3), k_E(4), k_E(5), k_E(6)], mu_s);    
    rr_pl(:, i) = r;
    EAA(i,1) =  acos(dot(y1(i,1:3),rr_pl(:,i)) / (norm(rr_pl(:,i))*norm(y1(i,1:3))));
    
end


SAA = zeros(length(t_IT),1);

for i = 1:length(t_IT)
       
    SAA(i,1) =  acos(dot(-y1(i,1:3),y1(i,4:6)) / (norm(y1(i,4:6))*norm(-y1(i,1:3)))) ;
    
end

figure(1)
plot(t_IT,rad2deg(EAA),'Linewidth', 2, 'DisplayName', 'Transfer arc', 'Color', 'r')
xlabel(' Time [day]')
ylabel('Angle [deg]')
title('Earth-Sun-s/c angle')
grid on; grid minor;

figure(2)
plot(t_IT,rad2deg(SAA),'Linewidth', 2, 'DisplayName', 'Transfer arc', 'Color', 'r')
xlabel(' Time [day]')
ylabel('Angle [deg]')
title('  Sun-s/c-\underline v angle ')
grid on; grid minor;

end

