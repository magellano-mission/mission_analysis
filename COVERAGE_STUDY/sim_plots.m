function sim_plots(SimType, wlk_vec, inclinations, bw, semi_major_axes, Min_cov_lat)

% bw, sma, incl, Nsat
% Colors:
% 1 = green = min 15 sats;
% 2 = yellow = min 18 sats;
% 3 = orange = min 21 sats;
% 4 = red = no coverage;

Ni = length(inclinations);
Na = length(semi_major_axes);
Nb = length(bw);

switch SimType
    case "Dani"
        
        C = (Min_cov_lat >= 4);
        for k = 1:Ni
            figure; hold on; grid on
            C_ith = squeeze(C(:, :, k, :));
            for i = 1:Nb
                for j = 1:Na
                    C_sat = C_ith(i, j, :);
                    if C_sat(1) == 1
                        h1 = plot(semi_major_axes(j), bw(i), 'og', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
                    elseif C_sat(1) == 0 && C_sat(2) == 1
                        h2 = plot(semi_major_axes(j), bw(i), 'oy', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
                    elseif C_sat(1) == 0 && C_sat(2) == 0 && C_sat(3) == 1
                        h3 = plot(semi_major_axes(j), bw(i), 'o', 'Color', [255, 191, 0]/255, 'MarkerSize', 10, 'MarkerFaceColor', [255, 191, 0]/255);
                    elseif C_sat(1) == 0 && C_sat(2) == 0 && C_sat(3) == 0
                        h4 = plot(semi_major_axes(j), bw(i), 'or', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
                    end
                end
            end
            
            if exist('h3', 'var') 
                legend([h1, h2, h3, h4], {'min 15 sat', 'min 18 sat', 'min 21 sat', 'no cov'}, 'Location', 'Best')
            else
               legend([h1, h2, h4], {'min 15 sat', 'min 18 sat', 'no cov'}, 'Location', 'Best') 
            end
            xlabel('SMA [km]'), ylabel('Beamwidth [deg]')
            title(strcat('Inclination = ',num2str(rad2deg(inclinations(k))), ' deg'))
        end

        
        
end


% figure1 = figure();
% sgtitle(strcat('Min coverage per latitude - f(SMA); N sats: ',num2str(walker(1)), ', Inc: ', num2str(rad2deg(inclin)),', N orbits: ',num2str(N_orbits)))
% LAT = linspace(lat(1), lat(2), disc(2));
%
% for j =1:length(bw)
% fix_bw = squeeze(Min_cov_lat(j,:,:));
% subplot(3,2,j); sss = pcolor(semi_major_axes, LAT, fix_bw' ); hold on
% sss.FaceColor = 'Interp';
% sss.EdgeColor = 'k';
% for cont = 1:length(LAT)
%     for cont2 = 1:length(semi_major_axes)
%         if fix_bw(cont2,cont)>=Treshold
%              plot(semi_major_axes(cont2), LAT(cont), 'r.','Markersize',10)
%         end
%     end
% end
%
% xlabel('SMA [km]'), ylabel('LAT [deg]')
% title(strcat('BW = ',num2str(bw(j)), ' deg')), colorbar
% end
%
% annotation(figure1,'textbox',...
%     [0.072875 0.925329428989751 0.14353125 0.0453879941434846],...
%     'String',{strcat('red dots: N sats visible >= ', num2str(Treshold))},...
%     'FitBoxToText','off');

end