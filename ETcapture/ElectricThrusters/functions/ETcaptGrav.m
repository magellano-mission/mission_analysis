function aJ2_car = ETcaptGrav(t, R, data)

% Perturbation due to gravity
if data.switchers.grav
    [a_J2, a_J3, a_J4] = GravPerturb(t, R, data);
else
    a_J2 = [0; 0; 0];
    a_J3 = [0; 0; 0];
    a_J4 = [0; 0; 0];
end

aJ2_car = a_J2 + a_J3 + a_J4;

end