%% CONTINUOS THRUST
close all, clear all, clc
%GAUSS INTEGRATOR REQUIRES SOME MODIFICATIONS
parameters.isPerturbed = 0;
parameters.isEP = 0;                                       % 1:EP thrust, 0 Chem thust
parameters.t_BO = 3000000;                                 % burnout time (chemical thrust)
parameters.T = [100; 0; 0];                                % thrust [N] (constant for the moment)
parameters.Isp = 200;                                       % specific impulse [s]
parameters.M0 = 5000;                                      % Total Mass of the s/c [kg]
parameters.c_r = 0.5;
parameters.Across_sun = 10;                                % Cross area related to the sun [m^2]
parameters.isInterp = 0;
parameters.t0sym = date2mjd2000([2021, 1, 1, 0, 0, 0]);    
parameters.tmax = date2mjd2000([2021, 1, 10, 0, 0, 0]);
parameters.event = 0;
parameters.opt = odeset('RelTol',1e-10, 'AbsTol',1e-10) ; %'Event', @event_cont_thrust);
%parameters.delta_v_req = 

if parameters.isInterp == 1
        mu = astroConstants(4);
        
        X0 = zeros(7,1);
        X0(1:6) = uplanet(parameters.t0sym, 3);
        [X0(1:3), X0(4:6)] = kep2car2(X0(1:6),mu);
    
        [T, states, parout] = cart_cont_thrust_model(X0, parameters);
        R_E = zeros(3,length(T));
        R_M = R_E;
        kep_E = zeros(length(T),6);
        kep_M = kep_E;

        for jj = 1:length(T)
            kep_E(jj,:) = uplanet(T(jj)/86400,3);
        end
        for jj = 1:length(T)
            rr = kep2car2(kep_E(jj,1:6), mu);
            R_E(:,jj) = rr;
        end

        for jj = 1:length(T)
            kep_M(jj,:) = uplanet(T(jj)/86400,4);
        end
        for jj = 1:length(T)
            rr = kep2car2(kep_M(jj,1:6), mu);
            R_M(:,jj) = rr;
        end

        r_ = zeros(3,length(T));
        figure()
        for jj = 1:6
            subplot(3,3,jj), 
            if jj >= 3
                states(:,jj) = rad2deg(states(:,jj));
                kep_E(:,jj) = rad2deg(kep_E(:,jj));
            end
        plot(T, states(:,jj), '.'), hold on,
        plot(T, kep_E(:,jj), 'r.'), hold on,    
        title(num2str(jj))
        end 
        subplot(3,3,[7 8 9]), plot(T,states(:,7)); title('propellant mass')

        figure()
         I = imread('Sun.jpg');                            % Mars image
        RI = imref2d(size(I));
        RI.XWorldLimits = [-180 180];                       % Mars image x sizes
        RI.YWorldLimits = [-90 90];                         % Mars image y sizes
        rM = almanac('Sun','Radius','kilometers','sphere');
        [X, Y, Z] = ellipsoid(0, 0, 0, rM, rM, rM, 100); % spheric centered Mars
        planet = surf(X, Y, -Z,'Edgecolor', 'none');
        hold on
        set(planet,'FaceColor','texturemap','Cdata',I)
        axis equal
        plot3(R_E(1,:), R_E(2,:), R_E(3,:)), hold on
        plot3(R_M(1,:), R_M(2,:), R_M(3,:),'r'), hold on

        for jj = 1:length(T)
            states(jj,3:6) = deg2rad(states(jj,3:6));
            rr = kep2car2(states(jj,1:6), mu);
            r_(:,jj) = rr;
        end
        kep0 = states(1,1:6);
        r_0 = zeros(3,360);
        for theta = 1:360
            kep0(6) = deg2rad(theta);
            rr0 = kep2car2(kep0, mu);
            r_0(:,theta) = rr0;
        end
        plot3(r_0(1,:),r_0(2,:),r_0(3,:),'y:'), hold on
        plot3(r_(1,:), r_(2,:), r_(3,:),'y'), hold on
        plot3(r_(1,1), r_(2,1), r_(3,1),'yo'), hold on
        set(gcf, 'color', 'k')
        set(gca, 'color', 'k','visible','off')

        figure()
        plot(T,states(:,7)); title('s/c mass')
else
        mu = astroConstants(14);
        J2 = 0.00196; %Mars J2
        X0 = zeros(7,1);
        X0(1:6) = [11500 0.1 deg2rad(25) 0.1 0.1 0]; 
        
        T_orb = 2 * pi * sqrt(X0(1)^3 / mu); % orbital period
        
        [T, states, parout] = gauss_cont_thrust_model(X0, parameters);
        r_ = zeros(3,length(T));
        
        figure()
        for jj = 1:6
            subplot(3,3,jj), 
            if jj >= 3
                states(:,jj) = rad2deg(states(:,jj));
            end
        plot(T, states(:,jj), '.'), hold on,
        title(num2str(jj))
        end 
        subplot(3,3,[7 8 9]), plot(T,states(:,7)); title('s/c mass')
        
        figure()
        
        I = imread('Mars.jpg');                            % Mars image
        RI = imref2d(size(I));
        RI.XWorldLimits = [-180 180];                       % Mars image x sizes
        RI.YWorldLimits = [-90 90];                         % Mars image y sizes
        rM = almanac('Mars','Radius','kilometers','sphere');
        [X, Y, Z] = ellipsoid(0, 0, 0, rM, rM, rM, 100); % spheric centered Mars
        planet = surf(X, Y, -Z,'Edgecolor', 'none');
        hold on
        set(planet,'FaceColor','texturemap','Cdata',I)
        axis equal

        
        for jj = 1:length(T)
            states(jj,3:6) = deg2rad(states(jj,3:6));
            rr = kep2car2(states(jj,1:6), mu);
            r_(:,jj) = rr;
        end
        kep0 = states(1,1:6);
        r_0 = zeros(3,360);
        for theta = 1:360
            kep0(6) = deg2rad(theta);
            rr0 = kep2car2(kep0, mu);
            r_0(:,theta) = rr0;
        end
        plot3(r_0(1,:),r_0(2,:),r_0(3,:),'y:'), hold on
        plot3(r_(1,:), r_(2,:), r_(3,:),'y'), hold on
        plot3(r_(1,1), r_(2,1), r_(3,1),'yo'), hold on
        set(gcf, 'color', 'k')
        set(gca, 'color', 'k','visible','off')
end
        



