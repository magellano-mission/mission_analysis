function [value, isterminal, direction] = ETCorrectionEvent(t, Y, data)
%{ 
EVENT: stop for the beginning of the maneuver
%}

stopCond = data.thetaStop;    % theta at which start the maneuver (from GA)
value(1) = mod(Y(6),2*pi) - stopCond;

if t > 1
    isterminal(1) = 1;
else 
    isterminal(1) = 0;
end

direction(1) = 0;