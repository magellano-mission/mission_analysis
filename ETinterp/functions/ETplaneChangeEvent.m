function [value, isterminal, direction] = ETplaneChangeEvent(~, Y, data)
%{ 
EVENT Function to stop the integration according to a change of inclination;
The inclination is computed here below and a negative inclination that has
to be setting in the config file is added.
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

h = cross(R, V);
i = acos(h(3)/norm(h));

value = i + data.i;
isterminal = 1;
direction = 0;
