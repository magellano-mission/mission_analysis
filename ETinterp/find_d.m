function errorTOF = find_d(d, ABC, muS, th_f, rMnorm, th_dotM, gammaM, TOF, n_int)
I = TOFcomp(d, ABC, muS, th_f, rMnorm, th_dotM, gammaM,n_int);
errorTOF = I - TOF*86400;
end