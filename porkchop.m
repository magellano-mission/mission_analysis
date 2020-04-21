function [minDV,Earth_time,Mars_time]=porkchop(in_date_min,fin_date_min,EM_SP,n_per)

[~, mu_s] = uplanet(0, 1);

in_date_max = date2mjd2000(in_date_min) +n_per*EM_SP/3600/24;
fin_date_max = date2mjd2000(fin_date_min) + n_per*EM_SP/3600/24;
 
t_in = date2mjd2000(in_date_min):EM_SP/50/3600/24:in_date_max;
t_fin = date2mjd2000(fin_date_min):EM_SP/50/3600/24:fin_date_max;


%% Lambert's problem: from Earth to Mars

DvEM = zeros(length(t_fin),length(t_in));

for i = 1 : length(t_in)
    
    tin = t_in(i);
    
    for j = 1 : length(t_fin)
        
        tfin = t_fin(j);
                     
        [k_E, ~] = uplanet(tin, 3);
        [r_E, v_E] = kep2car(k_E(1), k_E(2), k_E(3), k_E(4), k_E(5), k_E(6), mu_s);
        
        [k_M, ~] = uplanet(tfin, 4);
        [r_M, v_M] = kep2car(k_M(1), k_M(2), k_M(3), k_M(4), k_M(5), k_M(6), mu_s);
                
        [~ , ~ ,~, ~, VI, VF, ~, ] = lambertMR(r_E, r_M, (tfin- tin)*24*3600 , mu_s);
        
       Dv1 = norm(v_E - VI');
       Dv2 = norm(VF' - v_M);
        
        DvEM(j, i) = Dv2 + Dv1;        
        
    end
    
end

[minDV,pos_idx_min]=min(DvEM(:));

[row_idx_min,col_idx_min]=ind2sub(size(DvEM),pos_idx_min);

Earth_time = t_in(col_idx_min);
Mars_time = t_fin(row_idx_min);

%% Plot 
in_date_max = mjd20002date(in_date_max);

x_min = datenum(strcat(num2str(in_date_min(2)), '-', num2str(in_date_min(3)), '-', num2str(in_date_min(1))));
x_max = datenum(strcat(num2str(in_date_max(2)), '-', num2str(in_date_max(3)), '-', num2str(in_date_max(1))));
x = linspace(x_min, x_max, length(t_in));

fin_date_max = mjd20002date(fin_date_max);

y_min = datenum(strcat(num2str(fin_date_min(2)), '-', num2str(fin_date_min(3)), '-', num2str(fin_date_min(1))));
y_max = datenum(strcat(num2str(fin_date_max(2)), '-', num2str(fin_date_max(3)), '-', num2str(fin_date_max(1))));
y = linspace(y_min, y_max, length(t_fin));

figure(1)

contour(x, y, DvEM, 1:15)
title(' Preliminary Porkchop Plot ')
grid on; grid minor
colorbar

ax = gca;

ax.XAxis.TickValues = x_min : EM_SP/3600/24 : x_max;
ax.YAxis.TickValues = y_min : EM_SP/3600/24 : y_max;

xtickangle(45)
xlabel('Departure from Earth date')
xlim([x_min x_max])
datetick('x', 'mmm-yyyy', 'keepticks')

ytickangle(0)
ylabel('Arrival to Mars date')
ylim([y_min y_max])
datetick('y', 'mmm-yyyy', 'keepticks')

end
