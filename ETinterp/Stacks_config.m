%% CONFIGURATION FILE

% Specific Impulse (T in N)
fun_imp = @(X) 0.051346711071873*(X.*1e3).^2-6.382637268333952*(X.*1e3)+3554.166890279770;

% Power profile  (T in N)
fun_power = @(X) -0.039232229759447*(X.*1e3).^2 + 33.399284889333220*(X.*1e3) - 139.3692059725397;

% NS stacks
% dataNS1.Isp = 4300;                                       % specific impulse [s]
dataNS1.Isp = fun_imp;                                       % specific impulse [s]
dataNS1.Mdry = 1600;                                      % Total Mass of the s/c [kg]
dataNS1.n_int = 1000;
dataNS1.Tmax = 0.25; % N
dataNS1.n_int = 1000;

% NS stacks
% dataNS2.Isp = 4300;                                       % specific impulse [s]
dataNS2.Isp = fun_imp;                                       % specific impulse [s]
dataNS2.Mdry = 1600;                                      % Total Mass of the s/c [kg]
dataNS2.n_int = 1000;
dataNS2.Tmax = 0.25; % N
dataNS2.n_int = 1000;

% NS stacks
% dataNS3.Isp = 4300;                                       % specific impulse [s]
dataNS3.Isp = fun_imp;                                       % specific impulse [s]
dataNS3.Mdry = 1600;                                      % Total Mass of the s/c [kg]
dataNS3.n_int = 1000;
dataNS3.Tmax = 0.25; % N
dataNS3.n_int = 1000;

% RS1 stack
% dataRS1.Isp = 4300;                                       % specific impulse [s]
dataRS1.Isp = fun_imp;                                       % specific impulse [s]
dataRS1.Mdry = 1100;                                      % Total Mass of the s/c [kg]
dataRS1.n_int = 1000;
dataRS1.Tmax = 0.25; % N
dataRS1.n_int = 1000;

% RS2 stack
% dataRS2.Isp = 4300;                                       % specific impulse [s]
dataRS2.Isp = fun_imp;                                       % specific impulse [s]
dataRS2.Mdry = 1100;                                      % Total Mass of the s/c [kg]
dataRS2.n_int = 1000;
dataRS2.Tmax = 0.25; % N
dataRS2.n_int = 1000;

% ECS stack
% dataECS.Isp = 4300;                                       % specific impulse [s]
dataECS.Isp = fun_imp;                                       % specific impulse [s]
dataECS.Mdry = 1300;                                      % Total Mass of the s/c [kg]
dataECS.n_int = 1000;
dataECS.Tmax = 0.25; % N
dataECS.n_int = 1000;