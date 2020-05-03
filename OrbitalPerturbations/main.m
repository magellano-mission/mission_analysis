%% MAIN
clc, clear, 
close all

% Figure Initialization
set(0,'DefaultFigureUnits', 'normalized');
set(0,'DefaultFigurePosition',[0 0 1 1]);
set(0,'DefaultTextFontSize',18);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultAxesXGrid','on')
set(0,'DefaultAxesYGrid','on')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');


%%

H = [100:1000:25000];
% H = 8110.5; % altitude to have semimajor axis 11500 km
ACC_mag = zeros(length(H),9);
XX = [];

for j = 1:length(H)

h = H(j);

mu_m = astroConstants(14); %mu Mars
R_m = 3389.5; %radius Mars

% Orbital parametrs
a = R_m + h;
e = 0.00001;
i = deg2rad(55);
Om = 0;
om = 0;
theta = 0;

x_0 = [a, e, i, Om, om, theta]; %initial conditions

T = 2 * pi * sqrt(a^3 / mu_m); % orbital period

%t0 = date2mjd2000([2024,1,1,00,00,00]); %choose a date to start integration
t0 = 0;
t_span = [t0 t0+0.25*T]; % integration time

tic
[t, x] = ode113(@(t,x)gauss(t,x), t_span,  x_0);
toc

accelerations = zeros(length(t),9);

for kk = 1 : length(t)
    [~, accelerations(kk,:)] = gauss(t(kk), x(kk,:));
end

ACC_mag(j,:) = accelerations(1,:);

XX = [XX;x];  %keplerian coordinates from 0 to end

end

for k = 1:length(t) % correction to keep RAAN inside [-360 +360] range
    if (XX(k,4))>2*pi
        XX(k,4) = XX(k,4) - 2*pi;
    elseif (XX(k,4))<-2*pi
        XX(k,4) = XX(k,4) + 2*pi;
    end
end


figure
loglog(H, ACC_mag*1000, 'LineWidth', 1);
legend('Primary','J2','J3','J4','Phobos','Deimos','Sun','SRP')
title('Accelerations')
xlabel('altitude [km]')
ylabel('acceleration [m/s^2]')
grid on


Nrev = (60*60*24*365)/(T); %number of revolutions in one year 

%% Orbital parameters variation
% To perform it select only one altitude and simulate for a larger amount
% of time.
% Suggestion: set accelerations due to the moons (in gauss function) to 0 and do not integrate them
% otherwise it will never ends. Analytical ephemerides are needed to
% perform accurate long simulations.



% figure 
% subplot(1,3,1)
% plot(t/(60*60*24*30),(XX(:,1)))
% grid on
% xlabel('time [months]')
% xlim([0 120])
% ylabel('a [km]')
% subplot(1,3,2)
% plot(t/(60*60*24*30),rad2deg(XX(:,3)))
% grid on
% xlabel('time [months]')
% xlim([0 120])
% ylabel('inc [deg]')
% subplot(1,3,3)
% plot(t/(60*60*24*30),rad2deg(XX(:,4)))
% grid on
% xlabel('time [months]')
% xlim([0 120])
% ylabel('OM [deg]')

