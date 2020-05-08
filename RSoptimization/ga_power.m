
% UHF
% lb_sma = 3800;      ub_sma = 7500;
% lb_Guser = 0;       ub_Guser = 5;
% lb_bw = 50;         ub_bw = 70;

% X-band
lb_sma = 9000;       ub_sma = 12000;
lb_Guser = 10;       ub_Guser = 30;
lb_bw = 5;           ub_bw = 30;

lb = [lb_sma, lb_Guser, lb_bw];
ub = [ub_sma, ub_Guser, ub_bw];

[x, Ptot] = ga(@getSatellitePower, 3, [], [], [], [], lb, ub, @mycon);

%[~, Nsat] = getSatellitePower(x)