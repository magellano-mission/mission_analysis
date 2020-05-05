

lb_sma = 4800;      ub_sma = 4900;
lb_Guser = 0;       ub_Guser = 1;
lb_bw = 10;         ub_bw = 12;

lb = [lb_sma, lb_Guser, lb_bw];
ub = [ub_sma, ub_Guser, ub_bw];

x = ga(@getSatellitePower, 3, [], [], [], [], lb, ub, @mycon);