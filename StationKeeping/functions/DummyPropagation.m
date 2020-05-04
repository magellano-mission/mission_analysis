function [Cov_Results, time_map, Y] = DummyPropagation(data)

[Y, theta, ~, h_user] = const_orbits(data);
time_map = time_mapping(Y, theta, h_user, data);
[Cov_Results] = getMinCoverage(time_map);

end