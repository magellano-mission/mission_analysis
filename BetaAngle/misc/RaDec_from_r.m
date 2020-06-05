function [ra, dec] = RaDec_from_r(r)

l = r(1)/norm(r);
m = r(2)/norm(r);
n = r(3)/norm(r);

dec = asin(n);

if m > 0
    ra = acos(l/cos(dec));
else
    ra = 2*pi - acos(l/cos(dec));
end

end