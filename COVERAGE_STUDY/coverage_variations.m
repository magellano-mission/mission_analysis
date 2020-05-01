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

bw = 40;                    % NS beamwidth [deg]
% bw = 20;                  % RS beamwidth [deg]

% Parameters
lon = [-180, 180];          
lat = [-90,90];
disc = [10 10];             % discretization grid
alt = 0;                    % altitude at which evaluate the coverage ( ground level = 0)
timesteps = 1000;               
N_orbits = 3;               % number of orbits

%orbital periods computation
inclinations = 25:10:75;
semi_major_axes = 6500:1500:15500;
%%
Min_cov_lat = zeros(length(inclinations), length(semi_major_axes), disc(2));
for inc = 1:length(inclinations)
    INC = deg2rad(inclinations(inc));       % inclination [deg]
        for smax = 1:length(semi_major_axes)
            wbb = waitbar(((inc-1)*length(inclinations) + smax)/(length(inclinations)*length(semi_major_axes)));
            SMA = semi_major_axes(smax);              % semi-major axis [km]
            [YYY, T, THETA, H] = const_orbits(walker, bw, SMA, INC, timesteps, N_orbits, alt);
            [time_map, ~, ~] = time_mapping(walker, YYY, T, THETA, H, lon, lat, disc);
            [ ~, ~, cov] = getMinCoverage(disc, time_map);
            Min_cov_lat(inc, smax ,:) = cov'; 
        end
end


