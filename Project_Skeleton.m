%% Project Skeleton Code

clear; clc; close all;

%%initial parameter: unit: m, degree, rad/sec
r2 = 36/100; % cm  o2B
r3 = 60/100;
r6 = 120/100;
r8 = 8.4/100; 

%% Part 1- Calculations for kinematic variables, caculated based on loop closure eqn



%% Part 2 - Force and Moment Calculation


%initial parameters:
rrod = 0.005;

dtheta2 = 2 * 180/pi; % theta2 dot
ddtheta2 = 0; % theta2 doble-dot - second derivative

% Rod masses (calculated by hand)
m2 = 0.07634;
m3 = 0.2545;
m4 = 5;
m5 = 5;
m6 = 0.1272;

% Mass moments of intertia (calculated by hand)
Ig3 = 1/12*m3*(3*rrod^2 + r3^2);
Ig6 = 1/12*m6*(3*rrod^2 + 4*0.6^2);

% Define lists to be plotted
M12_list = [];
theta2_list = [];
F12_list = [];
F23_list = [];
F34_list = [];
F16_list = [];
F56_list = [];
F14_list = [];
F35_list = [];
Fs_list = [];
Ms_list = [];

% Initialize lists for plotting
theta2_list = 0:1:360;
theta3_list = zeros(size(theta2_list));
dtheta3_list = zeros(size(theta2_list));
ddtheta3_list = zeros(size(theta2_list));
theta6_list = zeros(size(theta2_list));
dtheta6_list = zeros(size(theta2_list));
ddtheta5_list = zeros(size(theta2_list));
r5_list = zeros(size(theta2_list));
dr5_list = zeros(size(theta2_list));
r7_list = zeros(size(theta2_list));
dr7_list = zeros(size(theta2_list));

for i = 1:length(theta2_list)
    theta2 = theta2_list(i);

    % Note my reference angles in the PDF I posted
    % Also, I changed my math to work with m instead of CM for CI units
    % Define angles and angular derivatives (derivatives done by hand)
    theta3 = 180 + asind((r2*sind(theta2) - r8) ./ r3);
    theta5 = theta3;
    theta6 = 180 - asind((r2*sind(theta2 + theta5) / r6)) - theta5;
    r5 = ((-r2 * sind(theta2 + theta6)) / sind(theta5 + theta6));
    dtheta3 = (r2*dtheta2*cosd(theta2))/(r3*cosd(theta3));
    dr5 = (r5*dtheta3*sind(theta3+theta6) + r2*dtheta2*sind(theta2+theta6)) / cosd(theta5+theta6);
    dtheta6 = (r5*dtheta3 + r2*dtheta2*cosd(theta2-theta5)) / (r6*cosd(theta6+theta5));
    ddtheta3 = (r2/r3)*(((-sind(theta2)*dtheta2^2 + cosd(theta2)*ddtheta2)*cosd(theta3)+cosd(theta2)*dtheta2*sind(theta3)*dtheta3)/cosd(theta3)^2);
    r7 = (r2 * cosd(theta2)) - (r3 * cosd(theta3));
    dr7 = (r2*dtheta2*cosd(theta2)) / (r3*cosd(theta3));

    dN = (dr5*dtheta3+r5*ddtheta3)+r2*(-sind(theta2-theta5)*(dtheta2-dtheta3)*dtheta2+cosd(theta2-theta5)*ddtheta2);
    N = r5*dtheta3 + r5*theta5*ddtheta3 + r2*(-sind(theta2-theta5)*(dtheta2-dtheta3)*dtheta2 + cosd(theta2-theta5)*ddtheta2);
    ddtheta6 = (cosd(theta6+theta5)*dN + sind(theta6+theta5)*(dtheta6+dtheta3)*N)/r6*cosd(theta6+theta5)^2;
  
    % Define Accelerations
    ag2x = -0.18*dtheta2^2*cosd(theta2);
    ag2y = -0.18*dtheta2^2*sind(theta2);
    ag3x = -0.36*dtheta2^2*cosd(theta2) - 0.60*(ddtheta3*sind(theta3) + dtheta3^2*cosd(theta3));
    ag3y = -0.36*dtheta2^2*sind(theta2) - 0.60*(ddtheta3*cosd(theta3) - dtheta3^2*sind(theta3));
    ag4 = -0.36*dtheta2^2*cosd(theta2) - 1.2*(ddtheta3*sind(theta3) + dtheta3^2*cosd(theta3));
    ag5x = -0.60*(ddtheta6*sind(theta6) + dtheta6^2*cosd(theta6));
    ag5y = 0.60*(ddtheta6*cosd(theta6) - dtheta6^2*sind(theta6));
    ag6x = -0.30*(ddtheta6*sind(theta6) + dtheta6^2*cosd(theta6));
    ag6y = 0.30*(ddtheta6*cosd(theta6) - dtheta6^2*sind(theta6));

    B = get_ma_vector( ...
        m2,m3,m4,m5,m6, ...
        ag2x,ag2y,ag3x,ag3y,ag4,ag5x,ag5y,ag6x,ag6y, ...
        ddtheta3, ddtheta6, ...
        Ig3, Ig6);
    A = get_A_matrix(theta2,theta3,theta6);

    x = A\ B; % Ax = B, solution for x; note that in MATLAB: A\B = B/A
    
    % Parse results
    F12x = x(1);
    F12y = x(2);
    F23x = x(3);
    F23y = x(4);
    F34x = x(5);
    F34y = x(6);
    F16x = x(7);
    F16y = x(8);
    F56x = x(9);
    F56y = x(10);
    F14 = x(11);
    F35 = x(12);
    M12 = x(13);

    % Shaking force and moment
    Fs = sqrt(F12x^2 + (F12y + F14)^2);
    Ms = M12 + (1.2*cosd(theta3) - 0.36*cosd(theta2))*F14;

    % Calculate magnitudes
    F12 = sqrt(F12x^2 + F12y^2);
    F23 = sqrt(F23x^2 + F23y^2);
    F34 = sqrt(F34x^2 + F34y^2);
    F16 = sqrt(F16x^2 + F16y^2);
    F56 = sqrt(F56x^2 + F56y^2);
    F14 = abs(F14);
    F35 = abs(F35);
    M12 = abs(M12);
    Ms = abs(Ms);

    % Add values to running list to plot later
    M12_list = [M12_list; M12];
    theta2_list(i) = theta2;
    F12_list = [F12_list; F12];
    F23_list = [F23_list; F23];
    F34_list = [F34_list; F34];
    F16_list = [F16_list; F16];
    F56_list = [F56_list; F56];
    F14_list = [F14_list; F14];
    F35_list = [F35_list; F35];
    Fs_list = [Fs_list; Fs];
    Ms_list = [Ms_list; Ms];

    % Store values for plotting
    theta3_list(i) = theta3;
    dtheta3_list(i) = dtheta3;
    ddtheta3_list(i) = ddtheta3;
    theta6_list(i) = theta6;
    dtheta6_list(i) = dtheta6;
    ddtheta5_list(i) = ddtheta3;
    r5_list(i) = r5;
    dr5_list(i) = dr5;
    r7_list(i) = r7;
    dr7_list(i) = dr7;

end

%% Plot Kinematic Variables
figure;
tiledlayout(5, 2);

titles = {'\theta_3 (deg)', 'd\theta_3 (deg/s)', 'dd\theta_3 (deg/s^2)', ...
          '\theta_6 (deg)', 'd\theta_6 (deg/s)', 'dd\theta_5 (deg/s^2)', ...
          'r_5 (m)', 'dr_5 (m/s)', 'r_7 (m)', 'dr_7 (m/s)'};

variables = {theta3_list, dtheta3_list, ddtheta3_list, theta6_list, dtheta6_list, ...
             ddtheta5_list, r5_list, dr5_list, r7_list, dr7_list};

for j = 1:length(variables)
    nexttile;
    plot(theta2_list, variables{j}, 'LineWidth', 1.5);
    grid on;
    title(titles{j}, 'Interpreter', 'tex');
    xlabel('\theta_2 (degrees)');
    ylabel(titles{j});
end

sgtitle('Kinematic Variables vs \theta_2');

% Regular and Polar plots:
% Might have to transpose the Force vectors for polar plot. Do so if needed
% Polar plot only works with radians so will have to do it accordingly

figure;
tiledlayout(5, 2); 

nexttile;
plot(theta2_list,M12_list)
grid on;
title('M_{12} vs \theta2')
xlabel('\theta_2   (degrees)')
ylabel('M_{12} (Nm)')

nexttile;
plot(theta2_list,F12_list)
grid on;
title('F_{12} vs \theta2')
xlabel('\theta_2   (degrees)')
ylabel('F_{12} (N)')

nexttile;
plot(theta2_list,F23_list)
grid on;
title('F_{23} vs \theta2')
xlabel('\theta_2   (degrees)')
ylabel('F_{23} (N)')

nexttile;
plot(theta2_list,F34_list)
grid on;
title('F_{34} vs \theta2')
xlabel('\theta_2   (degrees)')
ylabel('F_{34} (N)')

nexttile;
plot(theta2_list,F16_list)
grid on;
title('F_{16} vs \theta2')
xlabel('\theta_2   (degrees)')
ylabel('F_{16} (N)')

nexttile;
plot(theta2_list,F56_list)
grid on;
title('F_{56} vs \theta2')
xlabel('\theta_2   (degrees)')
ylabel('F_{56} (N)')

nexttile;
plot(theta2_list,F14_list)
grid on;
title('F_{14} vs \theta2')
xlabel('\theta_2   (degrees)')
ylabel('F_{14} (N)')

nexttile;
plot(theta2_list,F35_list)
grid on;
title('F_{35} vs \theta2')
xlabel('\theta_2   (degrees)')
ylabel('F_{35} (N)')

nexttile;
plot(theta2_list,Fs_list)
grid on;
title('F_{s} vs \theta2')
xlabel('\theta_2   (degrees)')
ylabel('F_{s} (N)')

nexttile;
plot(theta2_list,Ms_list)
grid on;
title('M_{s} vs \theta2')
xlabel('\theta_2   (degrees)')
ylabel('M_{s} (Nm)')
