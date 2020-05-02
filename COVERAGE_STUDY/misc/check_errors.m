function check_errors(SimType, wlk_vec, inclinations, bw, semi_major_axes)


if SimType == "FPCP" && (length(wlk_vec) ~= 1 || length(inclinations) ~= 1 || length(bw) ~= 1 || length(semi_major_axes) ~= 1)
    error('SimType FPCP not compatible with non-scalar parameters, check data.m');
end