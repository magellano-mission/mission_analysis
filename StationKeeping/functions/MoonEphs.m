function [Rmoon] = MoonEphs(mjd2000, data, Moon)

switch Moon
    case 'Phobos'
        Ephs = data.Ephs_Phobos;
        
    case 'Deimos'
        Ephs = data.Ephs_Deimos;
    
end

time = data.Phobos_time_mjd2000;

kep = zeros(6, 1);
for i = 1:6
    kep_table = Ephs(:, i);
    kep(i) = interp1(time, kep_table, mjd2000, 'linear');
end

Rmoon = kep2car_r_only(kep);
    

