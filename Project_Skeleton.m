%% Project Skeleton Code

clear; clc; close all;

%%initial parameter: unit: m, degree, rad/sec
r2 = 36/100; % cm  o2B
r3 = 60/100;
r6 = 120/100;
r8 = 8.4/100; 

%initial parameters:
rrod = 0.005;

dtheta2 = 2; % theta2 dot
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
theta2_list = 0:(2*pi/360):(2*pi);
theta3_list = zeros(size(theta2_list));
dtheta3_list = zeros(size(theta2_list));
ddtheta3_list = zeros(size(theta2_list));
theta6_list = zeros(size(theta2_list));
dtheta6_list = zeros(size(theta2_list));
ddtheta6_list = zeros(size(theta2_list));
r5_list = zeros(size(theta2_list));
dr5_list = zeros(size(theta2_list));
ddr5_list = zeros(size(theta2_list));
r7_list = zeros(size(theta2_list));
dr7_list = zeros(size(theta2_list));
ddr7_list = zeros(size(theta2_list));
a_list = zeros(size(theta2_list));

for i = 1:length(theta2_list)
    theta2 = theta2_list(i);

    % Define angles and angular derivatives (derivatives done by hand)
    theta3 = pi + asin((r2*sin(theta2) - r8) / r3);
    theta5 = theta3;
    theta6 = pi - asin((r2*sin(theta2 + theta5) / r6)) - theta5;
    dtheta3 = (r2*dtheta2*cos(theta2)) / (r3*cos(theta3));
    ddtheta3 = (r2/r3)*(((-sin(theta2)*dtheta2 + cos(theta2)*ddtheta2)*cos(theta3)+cos(theta2)*dtheta2*sin(theta3)*dtheta3)/cos(theta3)^2);

    % Define hand derived link position, velocity and acceleration equations
    r5 = (r2*cos(theta2) - r6*cos(theta6)) / cos(theta3); 
    dtheta6 = (-r5*dtheta3 + r2*dtheta2*cos(theta2-theta5)) / (r6*cos(theta6 + theta5));
    ddtheta6 = (-r2*sin(theta2+theta3)*(dtheta2+dtheta3)^2 + r2*cos(theta2+theta3)*(ddtheta2+ddtheta3) + r6*sin(theta6+theta3)*(dtheta6+dtheta3)^2) / (r6*cos(theta6+theta3)) - ddtheta3;
    dr5 = (r5*dtheta3*sin(theta3+theta6) + r2*dtheta2*sin(theta2+theta6)) / cos(theta5+theta6);
    ddr5 = (-r2*(dtheta2^2)*cos(theta2) + 2*dr5*dtheta3*sin(theta3) + r5*ddtheta3*sin(theta3) + r5*(dtheta3^2)*cos(theta3) + r6*ddtheta6*sin(theta6) + r6*(dtheta6^2)*cos(theta6)) / cos(theta3);
    r7 = (r2 * cos(theta2)) - (r3 * cos(theta3));
    dr7 = (r2*dtheta2*cos(theta2)) / (r3*cos(theta3));
    ddr7 = r3*(ddtheta3*sin(theta3) + (theta3^2)*cos(theta3)) - r2*(ddtheta2*sin(theta2) + (dtheta2^2)*cos(theta2));

    % Define Accelerations
    ag2x = -0.18*dtheta2^2*cos(theta2);
    ag2y = -0.18*dtheta2^2*sin(theta2);
    ag3x = -0.36*dtheta2^2*cos(theta2) - 0.60*(ddtheta3*sin(theta3) + dtheta3^2*cos(theta3));
    ag3y = -0.36*dtheta2^2*sin(theta2) - 0.60*(ddtheta3*cos(theta3) - dtheta3^2*sin(theta3));
    ag4 = -0.36*dtheta2^2*cos(theta2) - 1.2*(ddtheta3*sin(theta3) + dtheta3^2*cos(theta3));
    ag5x = -0.60*(ddtheta6*sin(theta6) + dtheta6^2*cos(theta6));
    ag5y = 0.60*(ddtheta6*cos(theta6) - dtheta6^2*sin(theta6));
    ag6x = -0.30*(ddtheta6*sin(theta6) + dtheta6^2*cos(theta6));
    ag6y = 0.30*(ddtheta6*cos(theta6) - dtheta6^2*sin(theta6));

    % Define variable for the Coriolis Effect
    a = abs(2*dr5*dtheta3);

    % Use defined function to get B of matrix
    B = get_ma_vector( ...
        m2,m3,m4,m5,m6, ...
        ag2x,ag2y,ag3x,ag3y,ag4,ag5x,ag5y,ag6x,ag6y, ...
        ddtheta3, ddtheta6, ...
        Ig3, Ig6);
    A = get_A_matrix(theta2,theta3,theta6);

    x = A \ B; % Ax = B, solution for x; note that in MATLAB: A\B = B/A
    
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
    Ms = M12 + (1.2*cos(theta3) - 0.36*cos(theta2))*F14;

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
    ddtheta6_list(i) = ddtheta6;
    r5_list(i) = r5;
    dr5_list(i) = dr5;
    ddr5_list(i) = ddr5;
    r7_list(i) = r7;
    dr7_list(i) = dr7;
    ddr7_list(i) = ddr7;
    a_list(i) = a;

end

%% Plot Kinematic Variables
figure;
tiledlayout(2, 2); % 2x2 grid for 4 combined plots

% --- Plot θ₃, dθ₃, ddθ₃ vs. θ₂ ---
nexttile;
plot(theta2_list, theta3_list, 'r', 'LineWidth', 1.5); hold on;
plot(theta2_list, dtheta3_list, 'g', 'LineWidth', 1.5);
plot(theta2_list, ddtheta3_list, 'b', 'LineWidth', 1.5);
grid on;
title('\theta_3, d\theta_3, dd\theta_3 vs \theta_2');
xlabel('\theta_2 (radians)');
ylabel('Value');
legend('\theta_3 (rad)', 'd\theta_3 (rad/s)', 'dd\theta_3 (rad/s^2)');
hold off;

% --- Plot θ₆, dθ₆, ddθ₆ vs. θ₂ ---
nexttile;
plot(theta2_list, theta6_list, 'r', 'LineWidth', 1.5); hold on;
plot(theta2_list, dtheta6_list, 'g', 'LineWidth', 1.5);
plot(theta2_list, ddtheta6_list, 'b', 'LineWidth', 1.5);
grid on;
title('\theta_6, d\theta_6, dd\theta_6 vs \theta_2');
xlabel('\theta_2 (radians)');
ylabel('Value');
legend('\theta_6 (rad)', 'd\theta_6 (rad/s)', 'dd\theta_6 (rad/s^2)');
hold off;

% --- Plot r₅ and dr₅ vs. θ₂ ---
nexttile;
plot(theta2_list, r5_list, 'r', 'LineWidth', 1.5); hold on;
plot(theta2_list, dr5_list, 'g', 'LineWidth', 1.5);
plot(theta2_list, ddr5_list, 'b', 'LineWidth', 1.5);
grid on;
title('r_5, dr_5, ddr_5 vs \theta_2');
xlabel('\theta_2 (radians)');
ylabel('Value');
legend('r_5 (m)', 'dr_5 (m/s)', 'ddr_5 (m/s^2)');
hold off;

% --- Plot r₇ and dr₇ vs. θ₂ ---
nexttile;
plot(theta2_list, r7_list, 'r', 'LineWidth', 1.5); hold on;
plot(theta2_list, dr7_list, 'g', 'LineWidth', 1.5);
plot(theta2_list, ddr7_list, 'b', 'LineWidth', 1.5);
grid on;
title('r_7, dr_7, ddr_7 vs \theta_2');
xlabel('\theta_2 (radians)');
ylabel('Value');
legend('r_7 (m)', 'dr_7 (m/s)', 'ddr_7 (m/s^s');
hold off;
sgtitle('Kinematic Variables vs \theta_2');


%% Plot Forces & Moments vs Input angle 
figure;
tiledlayout(5, 2); 

nexttile;
plot(theta2_list,M12_list)
grid on;
title('M_{12} vs \theta2')
xlabel('\theta_2   (radians)')
ylabel('M_{12} (Nm)')

nexttile;
plot(theta2_list,F12_list)
grid on;
title('F_{12} vs \theta2')
xlabel('\theta_2   (radians)')
ylabel('F_{12} (N)')

nexttile;
plot(theta2_list,F23_list)
grid on;
title('F_{23} vs \theta2')
xlabel('\theta_2   (radians)')
ylabel('F_{23} (N)')

nexttile;
plot(theta2_list,F34_list)
grid on;
title('F_{34} vs \theta2')
xlabel('\theta_2   (radians)')
ylabel('F_{34} (N)')

nexttile;
plot(theta2_list,F16_list)
grid on;
title('F_{16} vs \theta2')
xlabel('\theta_2   (radians)')
ylabel('F_{16} (N)')

nexttile;
plot(theta2_list,F56_list)
grid on;
title('F_{56} vs \theta2')
xlabel('\theta_2   (radians)')
ylabel('F_{56} (N)')

nexttile;
plot(theta2_list,F14_list)
grid on;
title('F_{14} vs \theta2')
xlabel('\theta_2   (radians)')
ylabel('F_{14} (N)')

nexttile;
plot(theta2_list,F35_list)
grid on;
title('F_{35} vs \theta2')
xlabel('\theta_2   (radians)')
ylabel('F_{35} (N)')

nexttile;
plot(theta2_list,Fs_list)
grid on;
title('F_{s} vs \theta2')
xlabel('\theta_2   (radians)')
ylabel('F_{s} (N)')

nexttile;
plot(theta2_list,Ms_list)
grid on;
title('M_{s} vs \theta2')
xlabel('\theta_2   (radians)')
ylabel('M_{s} (Nm)')

%% Plot corealis effect vs input angle 
figure;
plot(theta2_list,a_list)
grid on;
title('Coriolis Force vs \theta2')
xlabel('\theta_2 (radians)')
ylabel('Coriois Force (N)')
