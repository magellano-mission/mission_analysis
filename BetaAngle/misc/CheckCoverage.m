function [flag] = CheckCoverage(Cov_Results, data)

if isequal(data.lat, data.lat(Cov_Results.cov_lon >= data.trashold))
    flag = true;
else
    flag = false;
end
