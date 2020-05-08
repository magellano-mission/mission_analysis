function [lb, ub] = ETphasBoundaries(lb, ub, mi, r0, m, T, dth)

t = lb:100:ub;
N = length(t);
t_l = zeros(N, 1);

for i = 1:N
    t_l(i) = MarsETphasideal(t(i), r0, m, T, mi, dth);
end

plot(t, t_l);

if all(t_l > 0)
    error('no feasible points, try to raise ub')
else
    t_l = t_l(t_l > 0);
end

[~, ii] = min(t_l);

lb = t(ii)/2;
ub = t(ii)*3;


