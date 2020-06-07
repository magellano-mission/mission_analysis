function th_range = ThRange(R, V, data)

r = norm(R);
v = norm(V);
mi = data.mi;

% Computation of the eccentricity vector
E = ((v^2 - mi/r)*R - (dot(R, V)*V))/mi;
e = norm(E);

vr = dot(V, R)/r;
th = acos( max( min( dot(E, R)/(e*r) , 1) , -1) );

if vr < 0 && e > 1
    th = - th;
else 
    th = 2*pi - th;
end

if e + cos(th) > 0
    th_range = true;
else
    th_range = false;
end