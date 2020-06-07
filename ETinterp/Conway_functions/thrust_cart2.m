function [T_cart] = thrust_cart2( T, alpha, beta, th, r1vers, RCRRv)
% Transformation from reference plane to cartesian
rrplanvers = cos(th)*r1vers + sin(th)*cross(RCRRv, r1vers);
T_cart= T*cos(beta)*sin(alpha)*rrplanvers + T*cos(beta)*cos(alpha)*cross(RCRRv, rrplanvers) + T*sin(beta)*RCRRv;
end