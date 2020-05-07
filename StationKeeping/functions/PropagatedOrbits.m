function const = PropagatedOrbits(data)

% TT = data.walker(1);                            % total number of satellites
% Nplane = data.walker(2);                        % number of planes
% F = data.walker(3);                             % relative spacing between satellites in adjacent planes
% 
% n_sat = TT/Nplane;                              % satellites per plane

% Setup 
X0_kep = [data.sma, 1e-3, data.inc, 359.999/180*pi, pi/3, pi/3];      % initial condition

% Orbits computation
% Y = zeros(TT, data.NT, 6);                      % states matrix (3D)
% a = zeros(data.NT); e = a; i = a;
% RAAN = a; PA = a; th = a;
%  for ii = 1 : Nplane
%      for kk = 1 : n_sat
%          X0_kep(6) = (kk - 1) * pi * 2 /n_sat + (ii - 1) * F * 2 * pi / TT;
%          X0_kep(4) = (ii - 1) * 2 * pi / Nplane;
         [R, V]  = kep2car(X0_kep, data.mi); 
         X0 = [R; V];
         [T, Ysat] = ode113(@CarIntegration, [0, data.tend], X0, data.opt, data);
         for t = 1:length(Ysat(:, 1))
             [~, parout] = CarIntegration(T(t), Ysat(t, :), data);
             a(t) = parout{1};
             e(t) = parout{2};
             i(t) = parout{3};
             RAAN(t) = parout{4};
             PA(t) = parout{5};
             th(t) = parout{6};
         end
%          Y((ii-1)*n_sat + kk, :, :) = Ysat;
         
%      end
%  end
 const.T = T;
 const.Y = Ysat; const.a = a; const.e = e;  const.i = i;
 const.RAAN = RAAN;  const.PA = PA;  const.th = th;
end