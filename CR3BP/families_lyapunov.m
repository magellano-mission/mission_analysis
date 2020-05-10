    %% OF A PERIODIC ORBIT 
    %initialization
    clear, close all, clc

    mu = 0.0121506683;
    M1 =  -mu; 
    M2 = (1-mu);
    options = odeset('AbsTol', 1e-10, 'RelTol',1e-10);

    figure('name', 'Orbit')
%     plot3(M1, 0, 0, 'k.', 'markersize', 15), hold on
    plot3(M2, 0, 0, 'r.', 'markersize', 15), hold on,
    grid on, xlabel('x'), ylabel('y'), zlabel('z'), axis equal, view([0 0 90])
    
    step = (-0.01:0.005:0.08);
    %search of the periodic Lyapounov function
    X0_0 = [1.062; 0; 0; 0; 0.461; 0];
    X0 = X0_0;
    X01 = X0(1:3);
    T_half = 1.86;
    X0s = zeros(4,1);
    X0s(1:3) = X0(4:6);
    X0s(4) = T_half;
    opts = optimset('TolFun', 1e-13, 'TolX', 1e-13);
for dz = step
    wb = waitbar(dz/(step(end)-step(1)));
    X01(1) = X0_0(1) + dz; 
    [X0s, fval, ~, output] = fsolve(@(X) periodicity_search_Lyapunov(X, X01, mu, M1, M2, options), X0s, opts);
    X0new = [X01; X0s(1:3)];
    Tspan = [0 2*X0s(4)]; 
    [Tnew, ynew] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), Tspan, X0new, options);
    plot3(ynew(:,1), ynew(:,2), ynew(:,3), 'r'), hold on 
    plot3(ynew(1,1), ynew(1,2), ynew(1,3), 'k+'), hold on 
end
delete(wb)
%%  MONODROMY COMPUTATION
    clear, close all, clc

    mu = 0.0121506683;
    M1 =  -mu; 
    M2 = (1-mu);
    options = odeset('AbsTol', 1e-13, 'RelTol',1e-13);
    %search of the periodic Lyapounov function
    X0_0 = [1.062; 0; 0; 0; 0.461; 0];
    X0 = X0_0;
    X01 = X0(1:3);
    T_half = 1.86;
    X0s = zeros(4,1);
    X0s(1:3) = X0(4:6);
    X0s(4) = T_half;
    opts = optimset('TolFun', 1e-13, 'TolX', 1e-13,'Display','Off');
    [X0s, fval, ~, output] = fsolve(@(X) periodicity_search_Lyapunov(X, X01, mu, M1, M2, options), X0s, opts);
    X0new = [X01; X0s(1:3)];
    Tspan = [0 2*X0s(4)]; 
[T_p, y_p] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), Tspan, X0new, options);
%calcolo mondromy matrix
MDM = mondromoy_comp (1e-12, X0new, 2*X0s(4), 0, mu, M1, M2, options);
[V, D] = eig(MDM);
v = V(:,6);
display(diag(D))
figure()
for i = 1:length(y_p)
    phi_tau = mondromoy_comp (1e-12, X0new, T_p(i)+eps, 0, mu, M1, M2, options);
    X0 = y_p(i,:)' + 1e-5*phi_tau* v;
    [T_m, y_m] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), [0 10], X0, options);
    plot3(y_m(:,1), y_m(:,2), y_m(:,3), 'b'), hold on 
end
plot3(y_p(:,1), y_p(:,2), y_p(:,3), 'r'), hold on, axis equal
    
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
mu = 0.0121506683;
    M1 =  -mu; 
    M2 = (1-mu);
    options = odeset('AbsTol', 1e-13, 'RelTol',1e-13);
    %search of the periodic Lyapounov function
    X0_0 = [1.062; 0; 0; 0; 0.461; 0];
    X0 = X0_0;
    X01 = X0(1:3);
    T_half = 1.86;
    X0s = zeros(4,1);
    X0s(1:3) = X0(4:6);
    X0s(4) = T_half;
    opts = optimset('TolFun', 1e-13, 'TolX', 1e-13,'Display','Off');
    [X0s, fval, ~, output] = fsolve(@(X) periodicity_search_Lyapunov(X, X01, mu, M1, M2, options), X0s, opts);
    X0new = [X01; X0s(1:3)];
    Tspan = [0 2*X0s(4)]; 
    
    options = odeset('AbsTol', 1e-14, 'RelTol',1e-14);
    MONODROMY = reshape(y(end,7:42),[6,6]);
    
    [V,D] = eig(MONODROMY);
    EIGS = diag(D);
%     figure()
%     plot(real(EIGS), imag(EIGS), 'kx')
    
    cont = 1;
    figure()
    for jj = 6
        v = V(:,jj);
        for i=1:10:length(T)
        STM = reshape(y(i,7:42),[6,6]);
        Y0 = y(i,1:6)' + 1e-9*STM*v;
        [~, y_p] = ode45(@(t,x) CR3BP(t, x, mu, M1, M2), [0 10], Y0, options);
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
% %%
% figure ()
%     X0_0 = [1.062; 0; 0; 0; 0.461; 0];
%     X0 = X0_0;
%     X01 = X0(1:3);
%     T_half = 1.86;
%     X0s = zeros(4,1);
%     X0s(1:3) = X0(4:6);
%     X0s(4) = T_half;
%     stepp1 = linspace(step(1), step(end), 100);
%     dz = stepp1(2) - stepp1(1);
% for dz = stepp1
%     wb = waitbar(dz/(stepp1(end)-stepp1(1)));
%     X01(1) = X0_0(1) + dz; 
%     [X0s, fval, ~, output] = fsolve(@(X) periodicity_search_Lyapunov(X, X01, mu, M1, M2, options), X0s, opts);
%     X0new = [X01; X0s(1:3)];
%     Tspan = [0 2*X0s(4)]; 
%     % MONODROMY MATRIX COMPUTATION
%     phi = mondromoy_comp (1e-10, X0new, 2*X0s(4), 0, mu, M1, M2, options);
%     [V, D] = eig(phi);
%     DD = diag(D);
%     for jj = 1:6
%         if real(DD(jj)) < 0
%         X02 = X0new + 1e-6*V(:,jj);
%         [Tnew, ynew] = ode113(@(t,x) CR3BP(t, x, mu, M1, M2), 2*Tspan, X02, options);
%         coll = ['b','y', 'k','c','g','m'];
%         subplot(2,3,jj),plot3(ynew(:,1), ynew(:,2), ynew(:,3), coll(jj)), hold on 
%                         plot3(ynew(1,1), ynew(1,2), ynew(1,3), 'kx'), hold on, axis equal, view([0 0 90])
%         end
%     end
% end
% delete(wb)
%%
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

function dX = periodicity_search_Lyapunov(X0plan, X01, mu, E, M, options)
    X0 = [X01; X0plan(1:3)];
    Tspan = [0 X0plan(4)];
    [~, y] = ode113(@(t,X) CR3BP(t, X, mu, E, M), Tspan, X0, options);
    dX = [y(end,2); y(end,4); y(end,6); y(1,4); y(1,6)];
end

