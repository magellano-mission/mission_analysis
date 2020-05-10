    %% GENERATION OF A PERIODIC ORBIT 
    %initialization
    clear, close all, clc

    mu = 0.012150586550569;
    M1 =  -mu; 
    M2 = (1-mu);
    options = odeset('AbsTol', 1e-13, 'RelTol',1e-13);

    figure('name', 'Orbit')
%     plot3(M1, 0, 0, 'k.', 'markersize', 15), hold on
    plot3(M2, 0, 0, 'r.', 'markersize', 15), hold on,
    grid on, xlabel('x'), ylabel('y'), zlabel('z'), axis equal,
    
    stepp = -(0:0.005:0.15);
    %search of the periodic Lyapounov function
    X0_0 = [0.862; 0; 0.185; 0; 0.25; 0];
    X0 = X0_0;
    X03 = X0(3);
    T_half = 1.132;
    X0s = zeros(6,1);
    X0s(1:2) = X0(1:2);
    X0s(3:5) = X0(4:6);
    X0s(6) = T_half;
    opts = optimset('TolFun', 1e-13, 'TolX', 1e-13,'Display','Off');
for dz = stepp
    wb = waitbar(dz/(stepp(end)-stepp(1)));
    X03 = X0_0(3) + dz; 
    [X0s, fval, ~, output] = fsolve(@(X) periodicity_search_Halo(X, X03, mu, M1, M2, options), X0s, opts);
    X0new = [X0s(1:2); X03; X0s(3:5)];
    Tspan = [0 2*X0s(6)];
    [Tnew, ynew] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), Tspan, X0new, options);
    plot3(ynew(:,1), ynew(:,2), ynew(:,3), 'r'), hold on 
    plot3(ynew(1,1), ynew(1,2), ynew(1,3), 'k+'), hold on 
end
delete(wb)

%% MANIFOLD COMPUTATION
X0_0 = [0.862; 0; 0.185; 0; 0.25; 0];
X0 = X0_0;
X03 = X0(3);
T_half = 1.132;
X0s = zeros(6,1);
X0s(1:2) = X0(1:2);
X0s(3:5) = X0(4:6);
X0s(6) = T_half;
opts = optimset('TolFun', 1e-13, 'TolX', 1e-13,'Display','Off');
[X0s, fval, ~, output] = fsolve(@(X) periodicity_search_Halo(X, X03, mu, M1, M2, options), X0s, opts);
X0new = [X0s(1:2); X03; X0s(3:5)];
Tspan = [0 2*X0s(6)];
[T_p, y_p] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), Tspan, X0new, options);
%calcolo mondromy matrix
MDM = mondromoy_comp (1e-12, X0new, 2*X0s(6), 0, mu, M1, M2, options);
[V, D] = eig(MDM);
v = V(:,2);
display(diag(D))
figure()
plot3(y_p(:,1), y_p(:,2), y_p(:,3), 'r'), hold on 
plot3(y_p(1,1), y_p(1,2), y_p(1,3), 'k+'), hold on, axis equal
for i = 1:length(y_p)
    phi_tau = mondromoy_comp (1e-12, X0new, T_p(i)+eps, 0, mu, M1, M2, options);
    X0 = y_p(i,:)' + 1e-4 *phi_tau* v;
    [T_m, y_m] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), [0 10], X0, options);
    plot3(y_m(:,1), y_m(:,2), y_m(:,3), 'b'), hold on 
end

[V, D] = eig(MONODROMY);
v = V(:,2);
display(diag(D))
figure()
plot3(y_p(:,1), y_p(:,2), y_p(:,3), 'r'), hold on 
plot3(y_p(1,1), y_p(1,2), y_p(1,3), 'k+'), hold on, axis equal
    options = odeset('AbsTol', 1e-16, 'RelTol',1e-16);
for i = 1:length(y_p)
    phi_tau = mondromoy_comp (1e-12, X0new, T_p(i)+eps, 0, mu, M1, M2, options);
    X0 = y_p(i,:)' + 1e-4 *phi_tau* v;
    [T_m, y_m] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), [0 10], X0, options);
    plot3(y_m(:,1), y_m(:,2), y_m(:,3), 'b'), hold on 
end

%
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

   
    
    %% Manifold computation
    X0_0 = [0.862; 0; 0.185; 0; 0.25; 0];
    X0 = X0_0;
    X03 = X0(3);
    T_half = 1.132;
    X0s = zeros(6,1);
    X0s(1:2) = X0(1:2);
    X0s(3:5) = X0(4:6);
    X0s(6) = T_half;
    opts = optimset('TolFun', 1e-13, 'TolX', 1e-13,'Display','Off');
    [X0s, fval, ~, output] = fsolve(@(X) periodicity_search_Halo(X, X03, mu, M1, M2, options), X0s, opts);
    X0new = [X0s(1:2); X03; X0s(3:5)];
    Tspan = [0 2*X0s(6)];
    
    options = odeset('AbsTol', 1e-14, 'RelTol',1e-14);
    MONODROMY = reshape(y(end,7:42),[6,6]);
    
    [V,D] = eig(MONODROMY);
    EIGS = diag(D);
%     figure()
%     plot(real(EIGS), imag(EIGS), 'kx')
    
    cont = 1;
    figure()
    for jj = 4
        v = V(:,jj);
        for i=1:10:length(T)
        STM = reshape(y(i,7:42),[6,6]);
        Y0 = y(i,1:6)' + 1e-1*STM*v;
        [~, y_p] = ode45(@(t,x) CR3BP(t, x, mu, M1, M2), [0 2], Y0, options);
        colo = ['b','k','c','g'];
        plot3(y_p(:,1), y_p(:,2), y_p(:,3)), hold on,   
        [~, y_p0] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), Tspan, X0new, options);
        end
        plot3(y_p0(:,1), y_p0(:,2), y_p0(:,3),'r-', 'markersize',5), hold off, axis equal
        cont = cont + 1;
    end
    
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
% X0_0 = [0.862; 0; 0.185; 0; 0.25; 0];
%     X0 = X0_0;
%     X03 = X0(3);
%     T_half = 1.132;
%     X0s = zeros(6,1);
%     X0s(1:2) = X0(1:2);
%     X0s(3:5) = X0(4:6);
%     X0s(6) = T_half;
%     DD = zeros(6, length(stepp)); jj = 1;
% figure(), 
% for dz = stepp
%     wb = waitbar(dz/(stepp(end)-stepp(1)));
%     X03 = X0_0(3) + dz; 
%     [X0s, fval, ~, output] = fsolve(@(X) periodicity_search_Halo(X, X03, mu, M1, M2, options), X0s, opts);
%     X0new = [X0s(1:2); X03; X0s(3:5)];
%     Tspan = [0 2*X0s(6)]; 
%     % MONODROMY MATRIX COMPUTATION
%     phi = mondromoy_comp (1e-10, X0new, 2*X0s(6), 0, mu, M1, M2, options);
%     [V, D] = eig(phi);
%     DD(:,jj) = diag(D);
%     jj = jj +1;
% end
%     subplot(3,2,1), plot(real(DD(1,:)), imag(DD(1,:)), 'kx-'), hold on
%     subplot(3,2,2), plot(real(DD(2,:)), imag(DD(2,:)), 'kx-'), hold on
%     subplot(3,2,3), plot(real(DD(3,:)), imag(DD(3,:)), 'kx-'), hold on
%     subplot(3,2,4), plot(real(DD(4,:)), imag(DD(4,:)), 'kx-'), hold on
%     subplot(3,2,5), plot(real(DD(5,:)), imag(DD(5,:)), 'kx-'), hold on
%     subplot(3,2,6), plot(real(DD(6,:)), imag(DD(6,:)), 'kx-'), hold on
% delete(wb)

%
%     X0_0 = [0.862; 0; 0.185; 0; 0.25; 0];
%     X0 = X0_0;
%     X01 = X0(1);
%     T_half = 1.132;
%     X0s = zeros(6,1);
%     X0s(1) = X0_0(1);
%     X0s(2:5) = X0_0(3:6);
%     X0s(6) = T_half;
%     opts = optimset('TolFun', 1e-13, 'TolX', 1e-13);
% for dx = step
%     wb = waitbar(dx/(step(end)-step(1)));
%     X01 = X0_0(2) + dx; 
%     [X0s, fval, ~, output] = fsolve(@(X) periodicity_search_Halo1(X, X01, mu, M1, M2, options), X0s, opts);
%     X0new = [X0s(1); X01; X0s(2:5)];
%     Tspan = [0 2*X0s(6)]; 
%     [Tnew, ynew] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), Tspan, X0new, options);
%     plot3(ynew(:,1), ynew(:,2), ynew(:,3), 'g'), hold on 
%     plot3(ynew(1,1), ynew(1,2), ynew(1,3), 'k+'), hold on 
% end
% delete(wb)
function phi = mondromoy_comp (dx, X0, T_obs, T_phi,mu, M1, M2, options)
    Tspan = [T_phi T_obs];
    [~, y] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), Tspan, X0, options);
    phi = zeros(6,6);
    for i=1:6
        X0new = X0;
        X0new(i) = X0new(i) + dx;
        [~, ynew] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), Tspan, X0new, options);
        phi(:,i) = (ynew(end,:) - y(end,:))/dx;
    end
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
function dX = periodicity_search_Halo(X0plan, X03, mu, E, M, options)
    X0 = [X0plan(1:2); X03; X0plan(3:5)];
    Tspan = [0 X0plan(6)];
    [~, y] = ode113(@(t,X) CR3BP(t, X, mu, E, M), Tspan, X0, options);
    [~, y2] = ode113(@(t,X) CR3BP(t, X, mu, E, M), 2*Tspan, X0, options);
    dX = [y(end,2); y(end,4); y(end,6); y(1,4); y(1,6); y2(end,1)-y2(1,1); y2(end,2)-y2(1,2); y2(end,3)-y2(1,3)];
end
function dX = periodicity_search_Halo1(X0plan, X02, mu, E, M, options)
    X0 = [X0plan(1); X02; X0plan(2:5)];
    Tspan = [0 X0plan(6)];
    [~, y] = ode113(@(t,X) CR3BP(t, X, mu, E, M), Tspan, X0, options);
    [~, y2] = ode113(@(t,X) CR3BP(t, X, mu, E, M), 2*Tspan, X0, options);
    dX = [y(end,2); y(end,4); y(end,6); y(1,4); y(1,6); y2(end,1)-y2(1,1); y2(end,2)-y2(1,2); y2(end,3)-y2(1,3)];
end


