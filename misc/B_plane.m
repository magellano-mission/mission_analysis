function [BB, B, theta, STR] = B_plane(kep_hyp, rr,vv, mu)
%computation of B-plane
r = norm(rr); %position at soi
v = norm(vv); %v_inf

h_ = cross(rr,vv)/norm(cross(rr,vv));
ee = (v^2*rr - dot(r,v)*v)/mu - rr/r;
e = norm(ee);

beta = acos(1/e);

S_ = vv/v;

NN = [0;0;1];

T_ = cross(S_,NN)/norm(cross(S_,NN));

R_ = cross(S_,T_);

B = abs(kep_hyp(1))*sqrt(e^2 - 1);
B_ = cross(S_, h_);
BB = B*B_;

theta = acos(dot(BB,T_)/(norm(BB)*norm(T_)));

STR = [S_ T_ R_];
end

