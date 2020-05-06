function [RA, delta] = angles(r)
% r(j) : column vector
% RA = zeros(length(r),1);
RA = zeros(size(r,1),1);
rr = vecnorm(r,2,2);
l = r(:,1);
m = r(:,2);
n = r(:,3);
delta = asin(n./rr);
RA(m>0) = acos(l(m>0)./(rr(m>0).*cos(delta(m>0))));
% RA(m<0) = 2*pi- acos(l(m<0)./(rr(m<0).*cos(delta(m<0))));
RA(m<=0) = 2*pi- acos(l(m<=0)./(rr(m<=0).*cos(delta(m<=0))));
end