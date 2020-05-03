function sim_plots(SimType, wlk, inc, bw, sma, Min_cov_lat, lon, lat, N_min)

% bw, sma, incl, Nsat

Ni = length(inc);
Na = length(sma);
Nb = length(bw);

switch SimType
    case "Dani"
        System = ["NS", "RS"];
        Incl_string = string(round(inc*180/pi));
        
        C_ns = (Min_cov_lat >= 4);
        C_rs = (Min_cov_lat >= 1);
        for c = 1:2
            for k = 1:Ni
                fig = figure; hold on; grid on
                if c == 1
                    C = C_ns;
                elseif c == 2
                    C = C_rs;
                end
                
                C_ith = squeeze(C(:, :, k, :));
                for i = 1:Nb
                    for j = 1:Na
                        C_sat = C_ith(i, j, :);
                        if C_sat(1) == 1
                            h1 = plot(sma(j), bw(i), 'og', 'MarkerSize', 14, 'MarkerFaceColor', 'g');
                        elseif C_sat(1) == 0 && C_sat(2) == 1
                            h2 = plot(sma(j), bw(i), 'oy', 'MarkerSize', 14, 'MarkerFaceColor', 'y');
                        elseif C_sat(1) == 0 && C_sat(2) == 0 && C_sat(3) == 1
                            h3 = plot(sma(j), bw(i), 'o', 'Color', [255, 191, 0]/255, 'MarkerSize', 14, 'MarkerFaceColor', [255, 191, 0]/255);
                        elseif C_sat(1) == 0 && C_sat(2) == 0 && C_sat(3) == 0
                            h4 = plot(sma(j), bw(i), 'or', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
                        end
                    end
                end
                
                if exist('h3', 'var')
                    legend([h1, h2, h3, h4], {'min 15 sat', 'min 18 sat', 'min 21 sat', 'no cov'}, 'Location', 'Northeast')
                else
                    legend([h1, h2, h4], {'min 15 sat', 'min 18 sat', 'no cov'}, 'Location', 'Northeast')
                end
                Title = strcat(System(c), " coverage, i = ", Incl_string(k) );
                xlabel('sma [km]'), ylabel('beam-width [deg]'), title(Title)
                title_save = strcat(System(c), "_coverage_i=", Incl_string(k) );
                saveas(fig, strcat(title_save, '.png'))
                
            end
            
        end

    case 'FPCP'
        
        figure; grid off; hold on
        
        %mars texture reading
        texture = imread('Mars.jpg');
        x = [-180 180]; y = [-90 90];
        texture = flipud(texture);
        set(gca, 'YDir', 'reverse')
        
        colormap default
        image(x, y, texture), hold on,
        view(0,-90);
        axis([-180 180,-90,90])
        xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
        
        sss = pcolor(lon, lat, N_min'); hold on
        sss.FaceColor = 'Interp';
        sss.EdgeColor = 'interp';
        sss.FaceAlpha = 0.4;
        
        Rovers = {'Viking 1', 'Viking 2', 'Opportunity', 'Spirit', 'InSight', ...
            'Curiosity', 'Mars Polar Lander', 'ExoMars', 'Skipper (Mars-Penguin)'};
        Rov_lon = [-49.97, 134.29, -6, 175.47, 135.6, 137.4, 165.2, -24.3, -167];
        Rov_lat = [22.5, 47.6, -1.95, -14.57, 4.5, -4.6, -76.57, 18.1, -81];
        Text_lon = [3, 3, 3, -12, 3, 3, -40, 3, 3]; Text_lon = Text_lon + Rov_lon;
        Text_lat = [-3, -3, -3, -3, -3, -3, 5, -3, 3]; Text_lat = Text_lat + Rov_lat;
        
        for i = 1:length(Rovers)
            plot(Rov_lon(i), Rov_lat(i), 'wo', 'MarkerFace', 'r')        % Skipper (Mars-Penguin)
            text(Text_lon(i), Text_lat(i), Rovers(i), 'Color', 'w')
        end
        
        title('Minimum satellites coverage on surface', 'interpreter', 'latex', 'FontSize', 20)
        
        yline(0)
        xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
        ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
        colorbar
        
%     case 'a'
%         
%         figure1 = figure();
%         sgtitle(strcat('Min coverage per latitude - f(SMA); N sats: ',num2str(walker(1)), ', Inc: ', num2str(rad2deg(inclin)),', N orbits: ',num2str(N_orbits)))
%         LAT = linspace(lat(1), lat(2), disc(2));
%         
%         for j = 1:length(bw)
%             fix_bw = squeeze(Min_cov_lat(j,:,:));
%             subplot(3,2,j); sss = pcolor(sma, LAT, fix_bw' ); hold on
%             sss.FaceColor = 'Interp';
%             sss.EdgeColor = 'k';
%             for cont = 1:length(LAT)
%                 for cont2 = 1:length(sma)
%                     if fix_bw(cont2,cont)>=Treshold
%                         plot(sma(cont2), LAT(cont), 'r.','Markersize',10)
%                     end
%                 end
%             end
%             
%             xlabel('SMA [km]'), ylabel('LAT [deg]')
%             title(strcat('BW = ',num2str(bw(j)), ' deg')), colorbar
%         end
%         
%         annotation(figure1,'textbox',...
%             [0.072875 0.925329428989751 0.14353125 0.0453879941434846],...
%             'String',{strcat('red dots: N sats visible >= ', num2str(Treshold))},...
%             'FitBoxToText','off');

%% percent_cov, max_cov_gap, mean_cov_gap

% figure('Name', 'percent coverage')
% image(x, y, texture), hold on,
% view(0,-90);
% axis([-180 180,-90,90])
% xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
% sss1 = pcolor(lon, lat, percent_cov'); hold on
% sss1.FaceColor = 'Interp';
% sss1.EdgeColor = 'interp';
% sss1.FaceAlpha = 0.7;
% for i=1:length(lon)
%     for j=1:length(lat)
%         plot(lon(i),lat(j),'r+'), hold on
%     end
% end
% title('Percent coverage on surface', 'interpreter', 'latex', 'FontSize', 20)
% hold on
% yline(0)
% xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
% ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
% colorbar

%% 
% figure('Name','max coverage map')
% 
% image(x, y, texture), hold on,
% view(0,-90);
% axis([-180 180,-90,90])
% xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
% 
% sss = pcolor(lon, lat, max_cov_gap'); hold on
% sss.FaceColor = 'Interp';
% sss.EdgeColor = 'interp';
% sss.FaceAlpha = 0.7;
% 
% for i=1:length(lon)
%     for j=1:length(lat)
%         plot(lon(i),lat(j),'r+'), hold on
%     end
% end
% title('Maximum coverage gap', 'interpreter', 'latex', 'FontSize', 20)
% hold on
% yline(0)
% xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
% ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
% colorbar

%% 
% figure('Name','Mean coverage gap')
% 
% image(x, y, texture), hold on,
% view(0,-90);
% axis([-180 180,-90,90])
% xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
% 
% sss = pcolor(lon, lat, mean_cov_gap'); hold on
% sss.FaceColor = 'Interp';
% sss.EdgeColor = 'interp';
% sss.FaceAlpha = 0.7;
% 
% for i=1:length(lon)
%     for j=1:length(lat)
%         plot(lon(i),lat(j),'r+'), hold on
%     end
% end
% title('Minimum satellites coverage on surface with 40 deg beamwidth', 'interpreter', 'latex', 'FontSize', 20)
% hold on
% yline(0)
% xlabel('Longitude [deg]', 'interpreter', 'latex', 'FontSize', 15)
% ylabel('Latitude [deg]' , 'interpreter', 'latex', 'FontSize', 15)
% colorbar

%%
% figure
% axis equal
% 
% image(x, y, texture), hold on,
% view(0,-90);
% axis([-180 180,-90,90])
% xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
% timevarying = pcolor(lon, lat, time_map(:,:,1)');
% colorbar
% timevarying.FaceColor = 'Interp';
% timevarying.EdgeColor = 'interp';
% timevarying.FaceAlpha = 0.2;
% % timevarying = plot(lonSS(1), latSS(1),'x');
% 
% for jjjj = 1:timesteps
%     %     timevarying.XData = lonSS(jjjj); hold on
%     %     timevarying.YData = latSS(jjjj); hold on
%     %     xlim([lon(1) lon(2)])
%     %     ylim([lat(1) lat(2)])
%     timevarying.CData =  time_map(:,:,jjjj)';
%     drawnow
%     pause(0.07)
% end
% 
%% Constellation plot & groundtrack
% figure
% 
% mu = astroConstants(14);
% TT = walker(1);
% P = walker(2);
% F = walker(3);
% n_sat = TT/P;
% n_orbits = P;
% rM = almanac('mars','radius','kilometers','sphere');
% 
% subplot(1,2,1),
% I = imread('Mars.jpg');                            % Mars image
% RI = imref2d(size(I));
% RI.XWorldLimits = [-180 180];                       % Mars image x sizes
% RI.YWorldLimits = [-90 90];                         % Mars image y sizes
% [X, Y, Z] = ellipsoid(0, 0, 0, rM, rM, rM, 100); % spheric centered Mars
% planet = surf(X, Y, -Z,'Edgecolor', 'none');
% hold on
% set(planet,'FaceColor','texturemap','Cdata',I)
% axis equal
% set(gca,'Color','black')
% 
% subplot(1,2,2),
% 
% x=[-180 180]; y=[-90 90];
% texture = flipud(texture);
% image(x, y, texture), hold on,
% view(0,-90);
% axis([-180 180,-90,90])
% xlabel('Longitude [deg]'), ylabel('Latitude [deg]')
% pointY = animatedline('Color','y','Marker','o','MarkerSize',100);
%
% for jjjj=1:timesteps
%     for ii = 1 : n_orbits
%         colors = ['g','y','b'];
%         for kk = 1 : n_sat
%             Y = squeeze(YYY((ii-1)*n_sat + kk,:,:));
%             theta = THETA((ii-1)*n_sat + kk,:);
%             
%             subplot(1,2,1),plot3(Y(:,1),Y(:,2),Y(:,3), colors(ii)), hold on
%             addpoints(pointY, Y(jjjj,1),Y(jjjj,2),Y(jjjj,3)), hold on
%             
%             subplot(1,2,2),
%             [latSS, lonSS, ~] = groundTrack(Y, T);
%             sca = scircle1(latSS(jjjj), lonSS(jjjj), rad2deg(theta(1)));
%             geo = geoshow(sca(:,1),sca(:,2)); hold on
%             plot(lonSS,latSS,strcat(colors(ii),'.'),'MarkerSize',0.5)
%         end
%     end
%     pause(0.1)
%     subplot(1,2,1)
%     % clearpoints(pointY)
% end
    case 'none'

end

end