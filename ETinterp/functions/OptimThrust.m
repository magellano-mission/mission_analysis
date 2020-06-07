function aT_car = OptimThrust(R, V, T, data)


r = norm(R);
v = norm(V);
mi = data.mi;

% Computation of the eccentricity vector
E = ((v^2 - mi/r)*R - (dot(R, V)*V))/mi;
e = norm(E);

% Computation of the true anomaly
vr = dot(V, R)/r;
th = acos( max( min( dot(E, R)/(e*r) , 1) , -1) );

if vr < 0
    th = - th;
end
th;
H = cross(R, V);
h = norm(H);

s = (h + mi*r/h)*cos(th) + mi*e*r;
c = h*sin(th);
alpha = atan2(s, c);                            % angle from T to r

aT_rsw = T/data.M*1e-3*[cos(alpha); sin(alpha); 0];

vt = h/r;
vr = sqrt(v^2 - vt^2);
g = atan2(vt, vr);                              % fligth path angle

% rsw 2 tnh Matrix
G = [ cos(g)  -sin(g) 0
       sin(g) cos(g) 0 
       0         0      1 ];
   
% tnh 2 Car Matrix
t_hat = V/norm(V);
h_hat = cross(R, V)/norm(cross(R, V));
n_hat = cross(h_hat, t_hat);
A = [t_hat, n_hat, h_hat];
   
aT_tnh = G*aT_rsw;


aT_car = A*aT_tnh; 