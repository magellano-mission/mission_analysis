function errorTOF = find_d(d, ABC, muS, th_f, rM, th_dotM, gammaM, TOF, n_int)
I = TOFcomp(d, ABC, muS, th_f, rM, th_dotM, gammaM,n_int);
errorTOF = TOF*86400 - I;
end