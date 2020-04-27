clear,clc, close all


options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14 );
mu = 398600;

cd = 2.1; % [non dimensional]
arearatio = 0.02; % [m^2/kg]

k0 = [24263 0.7037 deg2rad(51.17) pi/2 pi/2 0.45];
[r,v] = kep2car(k0, mu);
state = [r; v];

thours= 24;
tspan = [0 thours*3600];


%% Integration of perturbed motion 
tic
[T1, Y] = ode113(@(t,y) cartesian(y, cd, arearatio), tspan, state, options);
toc

tic
[T2, K] = ode113(@(t,k) gauss(k, cd, arearatio), tspan, k0, options);
toc

R1 = Y(:, [1:3]);
R2 = zeros(length(T2), 3);

for i = 1:length(T2)
    [r, ~] = kep2car(K(i,:),mu);
    R2(i,:) = r';
end

KK =zeros(length(Y),6);
for i = 1:length(Y)
    KK(i,:) = car2kep(Y(i,1:3), Y(i,4:6), mu);
end

ths = K(:,6);
while any(ths > 2*pi)
    ths(ths > 2*pi) =  ths(ths > 2*pi) - 2*pi;
end

K(:,6) = ths;



%Costnat initial value for plots.
t = ones(1,length(T2));

a = k0(1)*t;
e = k0(2)*t;
i = k0(3)*t;
OM = k0(4)*t;
om = k0(5)*t;

%% ------ Plots ----

figure('Name', 'Orbit Propagation')
Terra3d
hold on
plot3(R1(:,1), R1(:,2), R1(:,3));
% plot3(R2(:,1), R2(:,2), R2(:,3));

figure('Name', 'Semimajor Axis')
plot(T2,K(:,1))
hold on
plot(T2,a)
plot(T1,KK(:,1))
ylim([22000 26000])

figure('Name', 'Eccentricity')
plot(T2,K(:,2))
hold on
plot(T2,e)
plot(T1,KK(:,2))
ylim([0.6 0.9 ])


figure('Name', 'Inclination')
plot(T2,K(:,3))
hold on
plot(T2,i)
plot(T1,KK(:,3))


figure('Name', 'Right Ascension of Ascending Note')
plot(T2,K(:,4))
hold on
plot(T2,OM)
plot(T1,KK(:,4))


figure('Name', 'Argument of Perigee')
plot(T2,K(:,5))
hold on
plot(T2,om)
plot(T1,KK(:,5))

% figure('Name', 'Anomaly')
% plot(T2,K(:,6))
% hold on
% % plot(T1,KK(:,6))


%% Time specifications:

N = size(T2,1);
T = tspan(2)/(N-1);
Fs = 1/T;
X = fftshift(fft(K(:,1)));
f = [-N/2:N/2-1]*Fs/N;

%%Plot the spectrum:
figure();
stem(f,abs(X)/N);
xlabel('Frequency (in hertz)');
title('Magnitude Response');

Xabs= abs(X)/N;
S = sort(Xabs,'descend');
I = find(Xabs == S(2))
f1 = f(I(2))
t1 = (1/f1)/(3600*24)