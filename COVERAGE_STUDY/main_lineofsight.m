% Coverage: lines of sight analysis

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

walker =[21,3,2];         % full constellation - original
% walker = [1,1,1];           % one satellite plot
% walker = [15,3,2];        % reduced constellation

% Parameters
SMA = 11500;                % semi-major axis [km]
INC = deg2rad(55);          % inclination [deg]
lon = [-180, 180];          
lat = [-90,90];
disc = [10 10];             % discretization grid
alt = 0;                    % altitude at which evaluate the coverage ( ground level = 0)
timesteps = 100000;               
N_orbits = 10;               % number of orbits
bw = 40;

%orbital periods computation
tic
[YYY, T, ~, ~] = const_orbits(walker, bw, SMA, INC, timesteps, N_orbits, alt);
toc
%%
TT = walker(1);                     % total number of satellites
P = walker(2);                      % number of planes
F = walker(3);                      % relative spacing between satellites in adjacent planes
    
n_sat = TT/P;                       % satellites per plane
n_orbits = P;

refplane = 1;                       %chooses the first satellite on the walker plane selected
Y = squeeze(YYY(1 + (refplane-1)*n_sat,:,:));

for k = 1:n_orbits
    figure ()
    sgtitle(strcat('plane', num2str(refplane),'- plane', num2str(k))), 
    if refplane == k
        n0 = 2;
    end

    for iii = n0:n_sat
        Y2 = squeeze(YYY((iii - 1)*n_orbits + k,:,:));
        [~, s_norm, s_dir] = lineofsight(Y, Y2);

        subplot(n_sat, 2, 2 * (iii-1) + 1)
        plot(T, s_norm, 'k')

        subplot(n_sat, 2, 2 * iii),
        plot(T, s_dir, 'k')
        
    end
    subplot(n_sat,2,1)
    title('norm variation')
    subplot(n_sat,2,2)
    title('direction wrt S0 variation')
    n0 = 1;
end

