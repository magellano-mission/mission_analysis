function check_errors(data)


if data.SimType == "FPCP" && (data.Nw ~= 1 || data.Ni ~= 1 || data.Nb ~= 1 || data.Na ~= 1)
    error('SimType FPCP not compatible with non-scalar parameters, check data.m');
end