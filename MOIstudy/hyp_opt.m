function [e_cap, dv_min, delta_opt] = hyp_opt(rp_desired, vv_inf, mu)

% {
% e_cap : [0 - 1] eccentricity of the capture orbit
% vv_inf : vector of the infinite velocity
% mu : gravity constant
% }

v_inf = norm(vv_inf);
e_cap = (2*mu/v_inf^2 - rp_desired)/(2*mu/v_inf^2 + rp_desired);
dv_min = v_inf * sqrt( (1-e_cap)/2 );
delta_opt = rp_desired * sqrt(2/(1-e_cap));

end

