function [t_box] = ReconstructTime(time_map, data)

ones_map = time_map == 0;

i = 1;
flag = true;
while flag && i < data.NT
    flag = not(any(any(ones_map(:, :, i))));
    i = i + 1;
end

t_box = ceil(data.tspan(i)/data.T_orb);