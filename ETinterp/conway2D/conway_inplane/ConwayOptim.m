function [y, m, T, r, TH, Vinf] = ConwayOptim(x, data_stacks)

t0 = x(1);              
TOF = x(2);
N_rev = x(3);
v_inf = x(4); 
alpha = x(5);

[kepE, miS] = uplanet(t0, 3);
[kepM, ~]    = uplanet(t0 + TOF, 4);

% rotating the initial velocity to get in Mars plane
[r1, v1] = kep2car2(kepE, miS); 
[r2, v2] = kep2car2(kepM, miS); 

%definition of plane of motion
r1norm = norm(r1);
r2norm = norm(r2);

r1vers = r1/r1norm;
r2vers = r2/r2norm;

hh = cross(r1vers, r2vers);
hvers = hh/norm(hh);
% dot(v1, hvers)
Vinf =  abs(dot(v2, hvers));
if Vinf > 0.1
    y = Inf;
else
    
    if hh(3) < 0
        hvers = -hvers;
    end
    v1 = v1 - dot(v1, hvers)*hvers;
    v2 = v2 - dot(v2, hvers)*hvers;
    
    v1 = v1 + v_inf*(cos(alpha)*r1vers + sin(alpha)*cross(hvers, r1vers));
    
    [T, m, r, TH] = ConwayInPlane(TOF, N_rev, r1norm, r2norm, r1vers, r2vers, hvers, v1, v2, miS, data_stacks);
    y = max(abs(T));
    
end


