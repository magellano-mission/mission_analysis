function [value, isterminal, direction] = ETsmaEvent(~, Y, data)
%{ 
EVENT Function to stop the integration according to a change of sma;
The sma is computed here below and subtracted by the starting sma and the 
needed change, both set in the config.
%}

% States
R = Y(1:3);

% input checks

value = norm(R) - data.sma - data.sma_change;
isterminal = 1;
direction = 0;
