function [box] = ReconstructBox(t_box, data)

box.a = data.adot*t_box;
box.e = data.edot*t_box;
box.i = data.idot*t_box*180/pi;