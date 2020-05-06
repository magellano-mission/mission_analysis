function check_errors(data)

S = data.SimType;
Nw = data.Nw;
Nb = data.Nb;
Na = data.Na;
Ni = data.Ni;

if data.Nsat/data.Norb ~= floor(data.Nsat/data.Norb)
    error('The number of satellites has to be a multiple of the number of walker planes, check data.m');
end

if S == "FPCP" && (Ni*Nw*Na*Nb ~= 1)
    error('SimType FPCP not compatible with non-scalar parameters, check data.m');
end

if S == "percent_cov" && (Ni*Nw*Na*Nb ~= 1)
     error('SimType percent_cov not compatible with non-scalar parameters, check data.m');
end

if S == "max_time_gap" && (Ni*Nw*Na*Nb ~= 1)
     error('SimType percent_cov not compatible with non-scalar parameters, check data.m');
end

if S == "mean_time_gap" && (Ni*Nw*Na*Nb ~= 1)
     error('SimType percent_cov not compatible with non-scalar parameters, check data.m');
end

if S == "plot_belli" && (Nw ~= 1 || Ni ~= 1)
    error('SimType plot_belli not compatible with non-scalar walker and inclinations, check data.m');
end



