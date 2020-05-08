function [Cov_Results] = getMinCoverage(time_map, data, tspan)

% disc: [discret_lon discret_lat]
% time_map: temporal-spatil history

n_lon = length(data.lon);
n_lat = length(data.lat);

% Initialization
percent_cov = zeros(n_lon, n_lat);
Cov_Results.max_cov_gap = zeros(n_lon, n_lat);
Cov_Results.mean_cov_gap = zeros(n_lon, n_lat);

% Computation of minimum
Cov_Results.N_min = min(time_map, [], 3);
Cov_Results.cov_lon = min(Cov_Results.N_min);
Cov_Results.cov_mars = min(Cov_Results.cov_lon);

% Computation of mean
Cov_Results.N_mean = mean(time_map, 3);



%percent coverage
for jjj = 1:data.NT
    for lo = 1 : n_lon
        for la = 1 : n_lat
            if time_map(lo, la, jjj) >= data.trashold
                percent_cov(lo, la) = percent_cov(lo, la) + 1;
            end
        end
    end
end
Cov_Results.percent_cov = percent_cov/data.NT*100;

cov_logic = (time_map >= data.trashold);
for lo = 1 : n_lon
    for la = 1 : n_lat
        a = squeeze(cov_logic(lo, la, :));
        t_cov = tspan(a)/3600;
        Ntcov = length(t_cov);
        if Ntcov >= 2
            gaps = t_cov(2:Ntcov) - t_cov(1:Ntcov-1);
            Cov_Results.max_cov_gap(lo, la) = max(gaps);
            Cov_Results.mean_cov_gap(lo, la) = mean(gaps);
        else
            Cov_Results.mean_cov_gap(lo, la) = tspan(end)/3600;
            Cov_Results.max_cov_gap(lo, la) = tspan(end)/3600;
        end
    end
end

end