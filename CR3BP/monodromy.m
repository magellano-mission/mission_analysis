
    %%
    options = odeset('AbsTol', 1e-13, 'RelTol',1e-13);
    tic
    [T, y] = ode113(@(t,x) myfun42(t,x, MM1, MM2), 2*Tspan, X00, options);
    toc
    figure()
    plot3(y(:,1), y(:,2), y(:,3)), axis equal
    %%
    function dx = myfun42(t,x, MM1, MM2)
    dx = zeros(42,1);
    dx(1:6) = MM2(t,x(1:6));
    dx(7:42)= MM1(t,x(7:42),x(1),x(2),x(3),x(4),x(5),x(6));
    end

function dx = CR3BP(t, x)

    d = sqrt((x(1) - E)^2 + x(2)^2 + x(3)^2);
    r = sqrt((x(1) - M)^2 + x(2)^2 + x(3)^2);
    %x [6,1] x, v
    dx = zeros(42,1);
    dx = zeros(6,1);
    dx(1:3) = x(4:6);
    dx(4) = x(1)+2*x(5)-(1-mu)*(x(1)+mu)/d^3-mu*(x(1)-1+mu)/r^3;
    dx(5) = x(2) - 2*x(4) - (1-mu)*x(2)/d^3 - mu*x(2)/r^3;
    dx(6) = -(1-mu)*x(3)/d^3 - mu*x(3)/r^3;
%     phi = zeros(6,6);
%     cont = 1;
%     for i = 1:6
%         for j = 1:6
%             phi(i,j) = x(cont);
%             cont = cont + 1;
%         end
%     end
%     A_tau = subs(A, xx, x(1));
%     A_tau = subs(A_tau, yy, x(2));
%     A_tau = subs(A_tau, zz, x(3));
%     A_tau = subs(A_tau, xd, x(4));
%     A_tau = subs(A_tau, yd, x(5));
%     A_tau = subs(A_tau, zd, x(6));
% 
%     dphi = A_tau*phi;
%     cont = 1;
%     for i = 1:6
%         for j = 1:6
%             dx(6 + cont) = dphi(i,j);
%             cont = cont + 1;
%         end
%     end
 x1 =@(x) x(1);
 x2 =@(x) x(2);
 x3 =@(x) x(3);
 dphi = MMM(t, x(7:42), x1, x2, x3);
 dx(7:42) = dphi;
    
end
function phi = wrapfun(x)
    phi = zeros(6,6);
    cont = 1;
    for i = 1:6
        for j = 1:6
            phi(i,j) = x(cont);
            cont = cont + 1;
        end
    end
end
function x = unwrapfun(phi)
    x = zeros(36,1);
    cont = 1;
    for i = 1:6
        for j = 1:6
            x(cont) = phi(i,j);
            cont = cont + 1;
        end
    end
end
% syms t x(t) y(t) z(t)
% xdot = diff(x,t);
% ydot = diff(y,t);
% zdot = diff(z,t);
% 
% xddot = diff(xdot,t);
% yddot = diff(ydot,t);
% zddot = diff(zdot,t);
% 
% %     autov = eig(A);
%     %%
%     AX = subs(A, xx, X0(1));
%     AX = subs(AX, yy, X0(2));
%     AX = subs(AX, zz, X0(3));
%     AX = subs(AX, xd, X0(4));
%     AX = subs(AX, yd, X0(5));
%     AX = subs(AX, zd, X0(6));
%     
%     phi0 = eye(6);
%     syms phi(t) [6,6]
%     dphi = diff(phi, t);

    