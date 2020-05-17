function distTH = MCI2TopocentricHorizon(distMCI, alpha, delta)

T = [-sin(alpha) cos(alpha) 0
     -sin(delta)*cos(alpha) -sin(delta)*sin(alpha) cos(delta)
     cos(delta)*cos(alpha) cos(delta)*sin(alpha) sin(delta)];


end