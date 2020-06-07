function [I, r, TH] = TOFcomp(d, ABC, muS, th_f, rMnorm, th_dotM, gammaM, n_int)
a = ABC(1); b = ABC(2); c = ABC(3);

A0 = [30*th_f^2   -10*th_f^3   th_f^4; ...
     -48*th_f    18*th_f^2  -2*th_f^3; ...
      20        -8*th_f      th_f^2];

A1 = [1/rMnorm-(a + b*th_f + c*th_f^2 + d*th_f^3); ...
    -tan(gammaM)/rMnorm-(b +2*c*th_f + 3*d*th_f^2); ...
     muS/(rMnorm^4*th_dotM^2)-(1/rMnorm + 2*c + 6*d*th_f)];
   
e = 0.5/th_f^6*A0(1,:)*A1;
f = 0.5/th_f^6*A0(2,:)*A1;
g = 0.5/th_f^6*A0(3,:)*A1;

TH = linspace(0, th_f, n_int);
TH2 = TH.^2; TH3 = TH.^3; TH4 = TH.^4; TH5 = TH.^5; TH6 = TH.^6;

r = 1./(a + b*TH + c*TH2 + d*TH3 + e*TH4 + f*TH5 + g*TH6); %[km]

%dtdtheta
dTOF  = ( abs(r.^4/muS .* (1./r + 2*c + 6*d*TH + 12*e*TH2 + 20*f*TH3 + 30*g*TH4)) ).^0.5; 

%dtheta
hh = TH(2) - TH(1);
%Cavalieri method
I = 0;
for i = 2:n_int
    tm = 0.5*(TH(i-1)+TH(i));
    rm = 1./(a + b*tm + c*tm^2 + d*tm^3 + e*tm^4 + f*tm^5 + g*tm^6);
    dTOFm = ( abs(rm^4/muS * (1/rm + 2*c + 6*d*tm + 12*e*tm^2 + 20*f*tm^3 + 30*g*tm^4)) )^0.5;
    I = I + hh/6*(dTOF(i-1) + 4*dTOFm + dTOF(i));
end
end