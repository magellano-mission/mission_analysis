function [dmdt] = massflowrate(time ,m, a, data_stacks)
T = a.*m*1000;
dmdt = -abs(T)/(9.81*data_stacks.Isp(T));
end