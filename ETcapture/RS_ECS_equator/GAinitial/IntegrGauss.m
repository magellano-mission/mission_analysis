function dY = IntegrGauss(~, Y, orbit)

a = orbit.a;
e = orbit.e;
T = orbit.T;
m0 = orbit.m0;
mi = orbit.mi;

i = Y(1);
th = Y(3);
om = Y(4);

ah = T/m0*1e-3;
p = a*(1 - e^2);
h = sqrt(mi*p);

r = p/(1 + e*cos(th));
thS = mod(th + om, 2*pi);

dY(1) = (r*cos(thS)/h)*ah;
dY(2) = r*sin(thS)*ah/(h*sin(i));
dY(3) = h/(r^2);
dY(4) = -r*sin(thS)*cos(i)/(h*sin(i))*ah;

dY = dY';
