% Coverage: Figures of Merit

%HP:
%circular orbits
%all constellation parameters fixed 
%
%timestep constant in all integrations
%groundtrack based approach
%conical shaped beam

%FoM: TO BE ADDED
% Percent coverage
% maximum gap
% mean gap
% time average gap
% mean response gap

% K_A = 2*pi;
% lambda0 = acos(rM/(norm(Y)));
% IAA = K_A*(1-cos(lambda0));


%% Setup
clear, close, clc

%walker = [18,3,2];         % full constellation - original
% walker = [1,1,1];           % one satellite plot
% walker = [15,3,2];        % reduced constellation

% bw = 40;                    % NS beamwidth [deg]
%bw = 20;                  % RS beamwidth [deg]
%bw = [20 40 70];
bw = 20;

% Parameters
lon = [-180, 180];          
lat = [-30,30];
disc = [15 6];             % discretization grid
alt = 0;                   % altitude at which evaluate the coverage ( ground level = 0)
timesteps = 1000;               
N_orbits = 3;              % number of orbits

% Ellipsoid model
para.rM_eq = 3393.4;         % equatorial radius [km]
rM_pol = 3375.7;        % polar radius [km]
para.ang_ecc = 2*atan(sqrt((para.rM_eq-rM_pol)/(para.rM_eq+rM_pol)));  % Mars angular eccentricity

T_lla2 = 1.02749125 * 24 * 3600;   % period [s]  
para.wM = 2 * pi / T_lla2;              % [rad/s]
para.theta_Airy_0 = 0;             % Mars principal meridian

%orbital periods computation
%inclinations = deg2rad(25:10:45);
inclinations = deg2rad(25);

%semi_major_axes = 6500:1000:11500;
semi_major_axes = 11500;

%wlk_vec = [15 18 21];
wlk_vec = 18;
Treshold = 4; 
%%
tic
%Min_cov_lat = zeros(length(bw), length(semi_major_axes), disc(2));
Min_cov_lat = zeros(length(bw), length(semi_major_axes),length(inclinations),length(wlk_vec));

for wlk_sat = 1:length(wlk_vec)
    walker = [wlk_vec(wlk_sat),3,2];
    for inc = 1:length(inclinations)
        INC = inclinations(inc);
        for beam = 1:length(bw)
            curr_beam = bw(beam);       % inclination [deg]
                for smax = 1:length(semi_major_axes)
                    %statebar = ((beam-1)*length(inclinations) + smax)/(length(inclinations)*length(semi_major_axes));
                    %waitforit = waitbar(statebar, strcat('Wait: ',num2str(statebar*100), '%'));
                    SMA = semi_major_axes(smax);              % semi-major axis [km]
                    [YYY, T, THETA, H] = const_orbits(walker, curr_beam, SMA, INC, timesteps, N_orbits, alt);
                    [time_map, ~, ~] = time_mapping(walker, YYY, T, THETA, H, lon, lat, disc, para);
                    [ ~, ~, cov] = getMinCoverage(disc, time_map);
                    cov = min(cov);     % minimo anche lungo latitudini
                    %Min_cov_lat(beam, smax ,:) = cov';
                    Min_cov_lat(beam, smax, inc, wlk_sat) = cov';
                    %delete(waitforit)
                end
        end
    end
end
toc

save('matriciozza.mat','Min_cov_lat')

return

% %% PLOT of performances
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
% 
