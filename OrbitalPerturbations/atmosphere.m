clear all
close all
clc
% Atmospheric model implemented in gauss. Data taken from MCD. 
% Just as a reference, no need to run.


% density data at Latitude 0.0N Longitude 0.0E Altitude 0.0 m ALS Local time 0.0h
% h : altitude[km]     r: density [kg/m^3]

h1 = 0;
r1 = 0.0150994;
h2 = 1;
r2 = 0.012456;
h3 = 5;
r3 = 0.00913664;
h4 = 10;
r4 = 0.00618213;
h5 = 20;
r5 = 0.00261632;
h6 = 30;
r6 = 0.00096723;
h7 = 40;
r7 = 0.000334501;
h8 = 50;
r8 = 0.000116851;
h9 = 60;
r9 = 4.1862e-05;
h10 = 70;
r10 = 1.40759e-05;
h11 = 80;
r11 = 4.073e-06;
h12 = 90;
r12 = 1.06172e-06;
h13 = 100;
r13 = 2.36307e-07;
h14 = 125;
r14 = 6.34858e-09;
h15 = 150;
r15 = 1.97726e-10;
h16 = 200;
r16 = 1.43765e-12;
h17 = 300;
r17 = 1.35052e-14;
h18 = 400;
r18 = 4.72894e-15;
h19 = 500;
r19 = 3.18361e-15;
h20 = 1000;
r20 = 7.69229e-16;



H1 = [h1,h2,h3,h4,h5,h6,h7,h8];
R1 = [r1,r2,r3,r4,r5,r6,r7,r8];

H2 = [h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19,h20];
R2 = [r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20];

coeff = fit(H1',R1','exp1');

a1 = 0.01442;
b1 = -0.08819;

coeff2 = fit(H2',R2','exp1');

a2 = 0.001039;
b2 = -0.04834;

density1 = @(x) a1*exp(b1*x);
density2 = @(x) a2*exp(b2*x);


altitude1 = linspace(h1,h8,1000);
altitude2 = linspace(h8,h20,1000);

figure
plot([H1,H2],[R1,R2],'o')
hold on
plot(altitude1,density1(altitude1));
plot(altitude2,density2(altitude2));


figure
plot(H2,R2,'o')
hold on
plot(altitude2,density2(altitude2));