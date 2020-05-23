function [rr, vv] = refplane2car( r, z,  vt, vr, th, r1vers, RCRRv)
% Transformation from reference plane to cartesian
rr = r*cos(th)*r1vers + r*sin(th)*cross(RCRRv, r1vers) + z*RCRRv;
vv = vr*r1vers + vt*cross(RCRRv, r1vers) + 0*RCRRv; %derivare
end