function f = elliptical_orbit(~, y, mu)
f = [y(4)
     y(5)
     y(6)
    -y(1)*mu / ( y(1)^2 + y(2)^2 + y(3)^2 )^(3/2) 
    -y(2)*mu / ( y(1)^2 + y(2)^2 + y(3)^2 )^(3/2)
    -y(3)*mu / ( y(1)^2 + y(2)^2 + y(3)^2 )^(3/2)];
end