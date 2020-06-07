function [T_i, alpha_i, beta_i, X_i, TOF_i] = interpolateDT(X, T, alpha, beta, TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, href, hh, v1, v2, muS, data)

[ ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, TOFdays] = ...
    Conway(TOF, N_rev, q, r1norm, r2norm, r1vers, r2vers, href, hh, v1, v2, muS, data);

T_i     = interp(T,50);
alpha_i = interp(alpha,50);
beta_i  = interp(beta,50);
TOF_i   = interp(TOFdays,50);
X_i     = zeros(3500,7);

for i = 1:7
    X_i(:,i) = interp(X(:,i),50);
end

figure(),
subplot(7,2,1), plot(X(:,1)), hold on, title('initial')
subplot(7,2,2), plot(X_i(:,1)),  title('interp')
subplot(7,2,3), plot(X(:,2)), 
subplot(7,2,4), plot(X_i(:,2)),  
subplot(7,2,5), plot(X(:,3)), 
subplot(7,2,6), plot(X_i(:,3)), 
subplot(7,2,7), plot(X(:,4)), 
subplot(7,2,8), plot(X_i(:,4)), 
subplot(7,2,9), plot(X(:,5)), 
subplot(7,2,10), plot(X_i(:,5)), 
subplot(7,2,11), plot(X(:,6)), 
subplot(7,2,12), plot(X_i(:,6)), 
subplot(7,2,13), plot(X(:,7)), 
subplot(7,2,14), plot(X_i(:,7)), 
figure()
subplot(3,2,1), plot(T), title('initial')
subplot(3,2,2), plot(T_i), title('interp')
subplot(3,2,3), plot(TOFdays, rad2deg(alpha)), title('$\alpha$'), ylabel('deg'), xlabel('TOF [days]')
subplot(3,2,4), plot(TOF_i, rad2deg(alpha_i)),title('$\alpha$'), ylabel('deg'), xlabel('TOF [days]')
subplot(3,2,5), plot(TOFdays, rad2deg(beta)),title('$\beta$'), ylabel('deg'), xlabel('TOF [days]')
subplot(3,2,6), plot(TOF_i, rad2deg(beta_i)),title('$\beta$'), ylabel('deg'), xlabel('TOF [days]')
end

