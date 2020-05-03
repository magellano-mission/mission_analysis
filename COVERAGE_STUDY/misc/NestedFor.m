function [Min_cov_lat, time_map, N_min, N_mean] = NestedFor(wlk_vec, inclinations, bw, semi_major_axes, lon, lat, timesteps, N_orbits, alt, para, perturb)

Min_cov_lat = zeros(length(bw), length(semi_major_axes),length(inclinations),length(wlk_vec));
for wlk_sat = 1:length(wlk_vec)
    walker = [wlk_vec(wlk_sat), 3, 2];
    for inc = 1:length(inclinations)
        INC = inclinations(inc);
        for beam = 1:length(bw)
            curr_beam = bw(beam);       % inclination [deg]
                for smax = 1:length(semi_major_axes)
                    SMA = semi_major_axes(smax);              % semi-major axis [km]
                    [YYY, T, THETA, ~, h_user] = const_orbits(walker, curr_beam, SMA, INC, timesteps, N_orbits, alt, perturb);
                    time_map = time_mapping(walker, YYY, T, THETA, h_user, lon, lat, para);
                    [N_min, N_mean, cov_lon, cov_mars] = getMinCoverage(lon, lat, time_map);
                    Min_cov_lat(beam, smax, inc, wlk_sat) = cov_mars;
                end
        end
    end
end