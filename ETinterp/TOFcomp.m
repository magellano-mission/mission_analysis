
function [I, r, TH] = TOFcomp(d, ABC, muS, th_f, rM, th_dotM, gammaM, n_int)
a = ABC(1); b = ABC(2); c = ABC(3);
rMnorm = norm(rM);
A0 = [30*th_f^2   -10*th_f^3   th_f^4; ...
     -48*th_f    18*th_f^2  -2*th_f^3; ...
      20        -8*th_f      th_f^2];

A1 = [1/rMnorm-(a + b*th_f + c*th_f^2 + d*th_f^3); ...
    -tan(gammaM)/rMnorm-(b +2*c*th_f + 3*d*th_f^2); ...
     muS/(rMnorm^4*th_dotM^2)-(1/rMnorm + 2*c + 6*d*th_f)];
   
e = 0.5/th_f^6*A0(1,:)*A1;
f = 0.5/th_f^6*A0(2,:)*A1;
g = 0.5/th_f^6*A0(3,:)*A1;

TH = linspace(0, th_f, n_int + 1);
for i = 1:n_int
    tm(i) = 0.5*(TH(i+1)+TH(i));
end

r = 1./(a + b*TH + c*TH.^2 + d*TH.^3 + e*TH.^4 + f*TH.^5 + g*TH.^6);
rm = 1./(a + b*tm + c*tm.^2 + d*tm.^3 + e*tm.^4 + f*tm.^5 + g*tm.^6);

dTOF = (abs(r.^4/muS .* (1./r + 2*c + 6*d*TH + 12*e*TH.^2 + 20*f*TH.^3 + 30*g*TH.^4)) ).^0.5;
dTOFm = (abs(rm.^4/muS .* (1./rm + 2*c + 6*d*tm + 12*e*tm.^2 + 20*f*tm.^3 + 30*g*tm.^4)) ).^0.5;

hh = TH(2) - TH(1);
%Cavalieri method
I = 0;
for i = 2:n_int+1
    I = I + dTOF(i-1) + 4*dTOFm(i-1) + dTOF(i);
end
I = I*hh/6;
end