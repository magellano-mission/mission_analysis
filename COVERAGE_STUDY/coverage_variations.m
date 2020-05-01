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

walker =[18,3,2];         % full constellation - original
% walker = [1,1,1];           % one satellite plot
% walker = [15,3,2];        % reduced constellation

% bw = 40;                    % NS beamwidth [deg]
bw = 20;                  % RS beamwidth [deg]

% Parameters
lon = [-180, 180];          
lat = [-90,90];
disc = [15 20];             % discretization grid
alt = 0;                    % altitude at which evaluate the coverage ( ground level = 0)
timesteps = 1000;               
N_orbits = 3;               % number of orbits

%orbital periods computation
inclinations = 25:10:75;
semi_major_axes = 6500:1000:15500;
Treshold = 4; 
%%
tic
Min_cov_lat = zeros(length(inclinations), length(semi_major_axes), disc(2));
for inc = 1:length(inclinations)
    INC = deg2rad(inclinations(inc));       % inclination [deg]
        for smax = 1:length(semi_major_axes)
            statebar = ((inc-1)*length(inclinations) + smax)/(length(inclinations)*length(semi_major_axes));
            waitforit = waitbar(statebar, strcat('Wait: ',num2str(statebar*100), '%'));
            SMA = semi_major_axes(smax);              % semi-major axis [km]
            [YYY, T, THETA, H] = const_orbits(walker, bw, SMA, INC, timesteps, N_orbits, alt);
            [time_map, ~, ~] = time_mapping(walker, YYY, T, THETA, H, lon, lat, disc);
            [ ~, ~, cov] = getMinCoverage(disc, time_map);
            Min_cov_lat(inc, smax ,:) = cov'; 
        end
end
toc
delete(waitforit)
%% PLOT of performances
figure1 = figure();
sgtitle(strcat('Min coverage per latitude - f(SMA); N sats:',num2str(walker(1)), ', Bw:', num2str(bw),', N orbits:',num2str(N_orbits)))
LAT = linspace(lat(1), lat(2), disc(2));
for j =1:length(inclinations)
fix_inc = squeeze(Min_cov_lat(j,:,:));
subplot(3,2,j); sss = pcolor(semi_major_axes, LAT, fix_inc' ); hold on
sss.FaceColor = 'Interp';
sss.EdgeColor = 'k';
for cont = 1:length(LAT)
    for cont2 = 1:length(semi_major_axes)
        if fix_inc(cont2,cont)>=1
             plot(semi_major_axes(cont2), LAT(cont), 'r.','Markersize',10)
        end
    end
end

xlabel('SMA [km]'), ylabel('LAT [deg]')
title(strcat('i = ',num2str(inclinations(j)), ' deg')), colorbar
end

annotation(figure1,'textbox',...
    [0.072875 0.925329428989751 0.14353125 0.0453879941434846],...
    'String',{strcat('red dots: N sats visible >= ', num2str(Treshold))},...
    'FitBoxToText','off');

