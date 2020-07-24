function [aT_car, thrust] = ETcaptThrust(t, R, V, M, data)

DayActual = data.InitDay + t/86400; 

KepS = uplanet(DayActual, 4);
rM2S = - kep2car_r_only(KepS);         % retriving the sun position vector wrt the planet
light = los(R, rM2S, data.R_pl);       % checking the line of sight

if light
    [aT_car, thrust] = ETConstThrust(rM2S, R, V, M, data);            % thrust [N]
    th_range = ThRange(R, V, data);

    if not(th_range)
        aT_car = [0; 0; 0];
        thrust = 0;
    end

else
    aT_car = [0; 0; 0];
    thrust = 0;
end



end