function [N_min, N_mean, cov] = getMinCoverage(LON, LAT, time_map)
%[N_min, N_mean, cov, percent_cov, max_cov_gap, mean_cov_gap] = getMinCoverage(disc, time_map, T_obs, trashold)

% disc: [discret_lon discret_lat]
% time_map: temporal-spatil history

n_lon = length(LON);
n_lat = length(LAT);

% Initialization
N_min = zeros(n_lon, n_lat);
N_mean = zeros(n_lon, n_lat);
%percent_cov = zeros(n_lon, n_lat);
%max_cov_gap = zeros(n_lon, n_lat);
% time_gap = zeros(n_lon, n_lat);
% time_gap_total = zeros(n_lon, n_lat);
% gap_count = zeros(n_lon, n_lat);
% 
% delta_t = T_obs/size(time_map,3);

% Computation of minimum
for lo = 1 : n_lon
    for la = 1 : n_lat
        N_min(lo,la) = min(time_map(lo,la,:)); % N_min = min(NN); % N_mesh(la, lo) = N_min;
    end
end

% Computation of mean
for lo = 1 : n_lon
    for la = 1 : n_lat
        N_mean(lo,la) = mean(time_map(lo,la,:));% N_min = min(NN); % N_mesh(la, lo) = N_min;
    end
end

N_mean = mean(time_map,3);

cov = min(N_min);

% %percent coverage
% for jjj = 1:size(time_map,3)
% for lo = 1 : n_lon
%     for la = 1 : n_lat
%         if time_map(lo,la,jjj)>=trashold
%             percent_cov(lo,la) = percent_cov(lo,la) + 1;
%         end
%     end
% end
% end
% percent_cov = percent_cov/size(time_map,3)*100;
% 
% %maximum time gap
% %first step
% jjj = 1;
% if time_map(lo,la,jjj)<trashold
%     time_gap(lo,la) = time_gap(lo,la) + 1;
%     time_gap_total(lo,la) = time_gap_total(lo,la) + 1;
% end
% if time_map(lo,la,jjj)>=trashold
%     time_gap(lo,la) = 0;
% end
% if time_gap(lo,la) > max_cov_gap(lo,la)
%     max_cov_gap(lo,la) = time_gap(lo,la);
% end
% 
% for jjj = 2:size(time_map,3)
% for lo = 1 : n_lon
%     for la = 1 : n_lat
%         if time_map(lo,la,jjj)<trashold
%             time_gap(lo,la) = time_gap(lo,la) + 1;
%             time_gap_total(lo,la) = time_gap_total(lo,la) + 1;
%         end
%         if time_map(lo,la,jjj)>=trashold
%             time_gap(lo,la) = 0;
%         end
%         if time_gap(lo,la) > max_cov_gap(lo,la)
%             max_cov_gap(lo,la) = time_gap(lo,la);
%         end
%         if (time_map(lo,la,jjj-1)>=trashold) && (time_map(lo,la,jjj)<=trashold)
%             gap_count(lo,la) = gap_count(lo,la) + 1;    
%         end
%     end
% end
% end
% 
% max_cov_gap = max_cov_gap*delta_t;
% mean_cov_gap = time_gap_total./gap_count*delta_t;

end