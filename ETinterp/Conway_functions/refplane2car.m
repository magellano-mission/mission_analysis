function [rr, vv] = refplane2car( r, z,  vt, vr, vz, th, r1vers, RCRRv)
% Transformation from reference plane to cartesian
rrplan = r*cos(th)*r1vers + r*sin(th)*cross(RCRRv, r1vers);
%the verical motion is just "added" on the planar motion, so i just took
%the radial versor in the plane and then added vz in h direction
% (it seems the right thing in Conway framework)
rr = rrplan + z*RCRRv;
rrplanvers = rrplan/(norm(rrplan));
vv = vr*rrplanvers + vt*cross(RCRRv, rrplanvers) + vz*RCRRv;
end