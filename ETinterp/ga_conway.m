function y = ga_conway(x, data_stacks)
t0 = x(1);              TOF = x(2);
N_rev = round(x(3));      q = x(4);

[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, m, T, ~] = ...
    Conway(t0, TOF, N_rev, q, data_stacks);

if ~isnan(m)

y(1) = max(abs(T));
y(2) = m(1) - m(end);

else
    y = 999999*(1*rand(2,1));
end
end

