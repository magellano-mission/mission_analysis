function errorTOF = find_d(d, ABC, muS, th_f, rM, th_dotM, gammaM, TOF)
I = TOFcomp(d, ABC, muS, th_f, rM, th_dotM, gammaM);
errorTOF = TOF*86400 - I;
end