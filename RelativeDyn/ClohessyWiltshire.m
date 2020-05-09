function[dX]=ClohessyWiltshire(t,X,par)
% Clohessy-Wiltshire : Relative dynamics- circular orbit - linearized
% approximation 
% X =[6x1]
% X(1:3) = r_rel 
% X(4:6) = v_rel 

% par    = struct: orbital parameters 

mu = par.mu; 
a  = par.a; % Circular orbit a = r 

% Initialize state vector 
dX = zeros(length(X),1); 
dX(1:3) = X(4:6);   % dr/dt 
n = sqrt(mu/a^3);   % orbital mean motion 

dX(4) = 3*n^2*X(1) +2*n*X(5);    %dv/dt
dX(5) = - 2*n*X(4);
dX(6) = -n^2*X(3);



