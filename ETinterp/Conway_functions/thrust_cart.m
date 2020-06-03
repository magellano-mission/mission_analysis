function [T_cart] = thrust_cart( T_in, T_out, gamma, th, r1vers, RCRRv)
% Transformation from reference plane to cartesian
rrplanvers = cos(th)*r1vers + sin(th)*cross(RCRRv, r1vers);
T_cart= T_in*sin(gamma)*rrplanvers + T_in*cos(gamma)*cross(RCRRv, rrplanvers) + T_out*RCRRv;
end