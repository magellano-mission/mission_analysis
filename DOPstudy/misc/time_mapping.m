function [time_map, GDOP_map, state] = time_mapping(curr, YYY, T, THETA, h_user, data)

% Setup
TT = curr.walker(1);                   % total number of satellites
Nplanes = curr.walker(2);              % number of planes
n_sat = TT/Nplanes;                    % saellites per plane

n_lon = length(data.lon);
n_lat = length(data.lat);

time_map = zeros(n_lon, n_lat, data.NT);
GDOP_map = zeros(n_lon, n_lat, data.NT);
state = cell(n_lon, n_lat, data.NT);
% HDOP_map = zeros(n_lon, n_lat, data.NT);
% VDOP_map = zeros(n_lon, n_lat, data.NT);

for jj = 1 : data.NT
    distM = zeros(TT, 3);
    curr_state = zeros(TT, 3);
    %state_cell = cell(n_lon, n_lat);
    instant_states = YYY(:, :, jj);
    instant_map = zeros(n_lon, n_lat);           % how many sats can view a certain point in a certain istant of time
    GDOP_inst_map = zeros(n_lon, n_lat);
%     HDOP_inst_map = zeros(n_lon, n_lat);
%     VDOP_inst_map = zeros(n_lon, n_lat);
    for lo = 1 : n_lon
        lon = data.lon(lo);                     % [deg]
        for la = 1 : n_lat          
            lat = data.lat(la);                 % [deg]
            % for each point, each satellite...
            for ii = 1 : Nplanes
                for kk = 1 : n_sat
                    %for each timestep, checking if point P is inside footprint
                    [kkk, rho, Yview] = coverageNumber(lon, lat, T(jj), instant_states((ii-1)*n_sat + kk, :),...
                        THETA((ii-1)*n_sat + kk, jj), h_user((ii-1)*n_sat + kk, jj), data);
                    instant_map(lo, la) = kkk + instant_map(lo, la);
                    if sum(rho) ~= 0 && sum(Yview) ~= 0
                        distM(instant_map(lo, la),:) = rho/norm(rho);
                        curr_state(instant_map(lo, la),:) = Yview;
                    end
                end
            end
            
            curr_state(all(curr_state == 0,2),:) = [];
            state{lo,la,jj} = curr_state;
            distM2 = [distM, -ones(TT,1)];
            distM2(all(distM2 == 0,2),:) = [];
            Q = inv(distM2'*distM2);
            GDOP_inst_map(lo,la) = sqrt(trace(Q));
%             HDOP_inst_map(lo,la) = sqrt(Q(1,1)+Q(2,2));
%             VDOP_inst_map(lo,la) = sqrt(Q(3,3));
        end
    end
    
    %state{:,:,jj} = state_cell;
    time_map(:,:,jj) = instant_map;
    GDOP_map(:,:,jj) = GDOP_inst_map;
    GDOP_map(time_map < 4) = 100;
%     HDOP_map(:,:,jj) = HDOP_inst_map;
%     HDOP_map(time_map < 4) = 100;
%     VDOP_map(:,:,jj) = VDOP_inst_map;
%     VDOP_map(time_map < 4) = 100;
end

end
