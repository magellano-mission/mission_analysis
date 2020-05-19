function dY = ThrustFreeInt(~, Y, data)
%{
ODE Function for the free dynamics
%}

% States
R = Y(1:3);
V = Y(4:6);

% input checks
if isrow(R)
    R = R';
end
if isrow(V)
    V = V';
end

% Derivative of the states
dY(1:3) = V;
dY(4:6) = - data.mi/norm(R)^3.*R;

dY = dY';

end