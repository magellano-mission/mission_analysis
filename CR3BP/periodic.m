%% GENERATION OF A PERIODIC ORBIT 
%initialization
clear, close all, clc

% %Halo Orbit
X0 = [0.862; 0; 0.185; 0; 0.25; 0];
T_half = 1.132;
Tspan = [0 2.5*T_half];
mu = 0.012150586550569;

% % Lyapunov Orbit
% X0 = [1.062; 0; 0; 0; 0.461; 0];
% T_half = 1.86;
% Tspan = [0 2*T_half];
% mu = 0.0121506683;

%%%%%%% 
%initial orbit
dist_ME = 1;
M1 = dist_ME *  -mu; 
M2 = dist_ME * (1-mu);
options = odeset('AbsTol', 1e-13, 'RelTol',1e-13);
[T, y] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), Tspan, X0, options);

%potential in the synodic rf 
xx = (0:.01:1.3);
yy = (-0.25:.01:0.25);
U = zeros(length(xx), length(yy));
for i=1:length(xx)
    for j = 1:length(yy)
        X = [xx(i); yy(j)];
        U(i,j) = U_computation(X, mu, M1, M2);
    end
end

% zero velocity curves 
syms xxx yyy
d = sqrt((xxx - M1)^2 + yyy^2); 
r = sqrt((xxx - M2)^2 + yyy^2); 

%guess orbit plot
figure('name', 'Orbit')
subplot(1,2,1), surf(xx, yy, U'/max(max(U)), 'facecolor', 'interp', 'edgecolor', 'interp'), hold on
axis equal 
subplot(1,2,2), plot3(y(:,1), y(:,2), y(:,3), 'k'), hold on 
plot3(M1, 0, 0, 'k.', 'markersize', 15), hold on
plot3(M2, 0, 0, 'r.', 'markersize', 15), hold on, 
JC = JC_computation(y, mu, M1, M2);
v2 = (xxx^2 + yyy^2) + 2*(1-mu)/d + 2*mu/r - mean(JC);
fcontour(v2, [-2 2 -2 2], 'LevelList', [0]), hold off
grid on, xlabel('x'), ylabel('y'), zlabel('z'), axis equal, title('initial guess')

% %search of the periodic orbit
% opts = optimset('TolFun', 1e-13, 'TolX', 1e-13);
% [X0s, fval, ~, output] = fsolve(@(X) periodicity_search(X, mu, M1, M2, options), [X0(1); X0(5); T_half], opts);
% display(fval)
% display(output.iterations)
% X0new = [X0s(1); 0; 0; 0; X0s(2); 0];
% Tspan = [0 2.5*X0s(3)]; 
% [Tnew, ynew] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), Tspan, X0new, options);
% JCnew = JC_computation(ynew, mu, M1, M2);
% figure()
% plot(Tnew, JCnew)

%search of the periodic orbit
opts = optimset('TolFun', 1e-13, 'TolX', 1e-13);
[X0s, fval, ~, output] = fsolve(@(X) periodicity_search(X, mu, M1, M2, options), [X0; T_half], opts);
display(fval)
display(output.iterations)
X0new = [X0s(1:6)];
Tspan = [0 5*X0s(7)]; 
[Tnew, ynew] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), Tspan, X0new, options);
JCnew = JC_computation(ynew, mu, M1, M2);
figure()
plot(Tnew, JCnew)

%plot of the new orbit
figure('name', 'Orbit')
plot3(ynew(:,1), ynew(:,2), ynew(:,3), 'k'), hold on 
plot3(M1, 0, 0, 'k.', 'markersize', 15), hold on
plot3(M2, 0, 0, 'r.', 'markersize', 15), hold on, 
grid on, xlabel('x'), ylabel('y'), zlabel('z'), view([0 0 90]), axis equal
v2 = (xxx^2 + yyy^2) + 2*(1-mu)/d + 2*mu/r - mean(JCnew);
fcontour(v2, [-2 2 -2 2], 'LevelList', [0]), hold off, title('periodic orbit')


%% computation of monodromy matrix

    syms xx yy zz xd yd zd xdd ydd zdd
    d = sqrt((xx - M1).^2 + yy.^2 + zz.^2);
    r = sqrt((xx - M2).^2 + yy.^2 + zz.^2);
    sys = [xd; ...
           yd; ...
           zd; ...
           xx + 2*yd-(1-mu)*(xx + mu)/d^3 - mu*(xx - 1+mu)/r^3; ...
           yy - 2*xd - (1-mu)*yy/d^3 - mu*yy/r^3; ...
           -(1-mu)*zz/d^3 - mu*zz/r^3];
       
    ss = [xx yy zz xd yd zd];
    A1 = diff(sys,xx);
    A2 = diff(sys,yy);
    A3 = diff(sys,zz);
    A4 = diff(sys,xd);
    A5 = diff(sys,yd);
    A6 = diff(sys,zd);
    A = [A1, A2, A3, A4, A5, A6];
    
    syms x(t) [6,1]
    A_tau = subs(A, xx, x1);
    A_tau = subs(A_tau, yy, x2);
    A_tau = subs(A_tau, zz, x3);
    A_tau = subs(A_tau, xd, diff(x1,t));
    A_tau = subs(A_tau, yd, diff(x2,t));
    A_tau = subs(A_tau, zd, diff(x3,t));
    
    sys = subs(sys, xx, x1);
    sys = subs(sys, yy, x2);
    sys = subs(sys, zz, x3);
    sys = subs(sys, xd, diff(x1,t));
    sys = subs(sys, yd, diff(x2,t));
    sys = subs(sys, zd, diff(x3,t));
    
    syms phi(t) [6,6]
    d = sqrt((x1 - M1).^2 + x2.^2 + x3.^2);
    r = sqrt((x1 - M2).^2 + x2.^2 + x3.^2);
    
    eqn6 = diff(x,t) == A_tau*x;
    dphi = A_tau*phi;
    dphi36 = reshape(diff(phi,t),[36,1]);
    Aphi36 = reshape(A_tau*phi,[36,1]);
    eqn8 = dphi36 == Aphi36;
    V1 = odeToVectorField(eqn8);
    V2 = odeToVectorField(eqn6);
    
    V1 = subs(V1, diff(x1,t), xd);
    V1 = subs(V1, diff(x2,t), yd);
    V1 = subs(V1, diff(x3,t), zd);
    V1 = subs(V1, x1, xx);
    V1 = subs(V1, x2, yy);
    V1 = subs(V1, x3, zz);
    
    V2 = subs(V2, diff(x1,t), xd);
    V2 = subs(V2, diff(x2,t), yd);
    V2 = subs(V2, diff(x3,t), zd);
    V2 = subs(V2, x1, xx);
    V2 = subs(V2, x2, yy);
    V2 = subs(V2, x3, zz);

    MM1 = matlabFunction(V1,'vars',{'t','Y','xx', 'yy', 'zz', 'xd', 'yd', 'zd'});
    MM2 = matlabFunction(V2,'vars',{'t','Y'});
    phi0 = reshape(eye(6),[36,1]);
    X00 = [X0new; phi0];
    options = odeset('AbsTol', 1e-13, 'RelTol',1e-13);
    [T, y] = ode113(@(t,x) myfun42(t,x, MM1, mu, M1, M2), Tspan, X00, options);
    figure()
    plot3(y(:,1), y(:,2), y(:,3)), axis equal
    
    MONODROMY = reshape(y(end,7:42),[6,6]);
    
    EIGS = eig(MONODROMY);
    figure()
    plot(real(EIGS), imag(EIGS), 'kx')
    %%
    function dx = myfun42(t,x, MM1, mu, E, M)
    d = sqrt((x(1) - E)^2 + x(2)^2 + x(3)^2);
    r = sqrt((x(1) - M)^2 + x(2)^2 + x(3)^2);
    dx = zeros(42,1);
    dx(1:3) = x(4:6);
    dx(4) = x(1)+2*x(5)-(1-mu)*(x(1)+mu)/d^3-mu*(x(1)-1+mu)/r^3;
    dx(5) = x(2) - 2*x(4) - (1-mu)*x(2)/d^3 - mu*x(2)/r^3;
    dx(6) = -(1-mu)*x(3)/d^3 - mu*x(3)/r^3;
    dx(7:42)= MM1(t,x(7:42),x(1),x(2),x(3),x(4),x(5),x(6));
    end

function dx = CR3BP(t, x, mu, E ,M)
    d = sqrt((x(1) - E)^2 + x(2)^2 + x(3)^2);
    r = sqrt((x(1) - M)^2 + x(2)^2 + x(3)^2);
    %x [6,1] x, v
    dx = zeros(6,1);
    dx(1:3) = x(4:6);
    dx(4) = x(1)+2*x(5)-(1-mu)*(x(1)+mu)/d^3-mu*(x(1)-1+mu)/r^3;
    dx(5) = x(2) - 2*x(4) - (1-mu)*x(2)/d^3 - mu*x(2)/r^3;
    dx(6) = -(1-mu)*x(3)/d^3 - mu*x(3)/r^3;
end
function dX = periodicity_search(X0plan, mu, E, M, options)
    X0 = [X0plan(1:6)];
    Tspan = [0 X0plan(7)];
    [~, y] = ode113(@(t,X) CR3BP(t, X, mu, E, M), Tspan, X0, options);
    [~, y2] = ode113(@(t,X) CR3BP(t, X, mu, E, M), 2*Tspan, X0, options);
    dX = [y(end,2); y(end,4); y(end,6); y(1,4); y(1,6); y2(end,1)-y2(1,1); y2(end,2)-y2(1,2); y2(end,3)-y2(1,3)];
end
% function dX = periodicity_search(X0plan, mu, E, M, options)
%     X0 = [X0plan(1); 0; X0plan(2); 0; X0plan(3); 0];
%     Tspan = [0 X0plan(4)];
%     [~, y] = ode113(@(t,X) CR3BP(t, X, mu, E, M), Tspan, X0, options);
%     dX = [y(end,2); y(end,4); y(end,6)];
%     
% end
function U = U_computation(x, mu, E ,M)
    x1 = x(1);
    x2 = x(2);
    d = sqrt((x1 - E).^2 + x2.^2); 
    r = sqrt((x1 - M).^2 + x2.^2);
    U =  0.5*(x1.^2+x2.^2) + (1-mu)/d + mu/r;
end
function JC = JC_computation(x, mu, E ,M)
    d = sqrt((x(:,1) - E).^2 + x(:,2).^2 + x(:,3).^2);
    r = sqrt((x(:,1) - M).^2 + x(:,2).^2 + x(:,3).^2);
    JC =  (x(:,1).^2+x(:,2).^2) + 2*(1-mu)./d + 2*mu./r - sum(abs(x(:,4:6).^2),2);
end

