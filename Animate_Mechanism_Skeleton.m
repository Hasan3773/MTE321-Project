%% Animate the motion of the mechanism 
clear; close all; clc;

% Define link lengths (units: cm)
r2 = 36; 
r3 = 120; 
r6 = 60; 
r8 = 8.4; 

% Define input angles (two full rotations)
inp_theta2 = 0:1:720; 

% Setup figure
figure(1);
hold on;
grid on;
axis([-75 200 -50 75]);
title('Mechanism Animation');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');

% Plot fixed pivot A
plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Origin');

% Initialize sliders and links to animate
link_AB = plot(nan, nan, 'ro-', 'LineWidth', 2, 'DisplayName', 'Crank (A → B)');
link_BD = plot(nan, nan, 'bo-', 'LineWidth', 2, 'DisplayName', 'Coupler (B → D)');
link_AC = plot(nan, nan, 'ko-', 'LineWidth', 2, 'DisplayName', 'Ground Link (A → C)');
slider_D = plot(nan, nan, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Slider D');
slider_C = plot(nan, nan, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8, 'DisplayName', 'Slider C');

legend('Location', 'northeastoutside');

% allow time for MATLAB to render
pause(0.5);

for i = 1:length(inp_theta2)

    % Define function for angles
    theta_2 = inp_theta2(i); 
    theta_3 = asind((r2 * sind(theta_2) - r8) / r3) + 180;
    theta_5 = theta_3;
    theta_6 = 180 - asind((r2 * sind(theta_2 + theta_5) / r6)) - theta_5;

    % Define point B as the crank
    Bx = r2 * cosd(theta_2);
    By = r2 * sind(theta_2);

    % Define hand derived functions for lengths
    r7 = (r2 * cosd(theta_2)) - (r3 * cosd(theta_3)); 
    r5 = ((-r2 * sind(theta_2 + theta_6)) / sind(theta_5 + theta_6));
    
    % Slider D as a function using lengths
    Dx = r7;
    Dy = r8;
    
    % Slider C using derived kinematic equations
    Cx = r6 * cosd(theta_6);
    Cy = r6 * sind(theta_6);

    % Update plots using new calculated values
    set(link_AB, 'XData', [0, Bx], 'YData', [0, By]);
    set(link_BD, 'XData', [Bx, Dx], 'YData', [By, Dy]);
    set(link_AC, 'XData', [0, Cx], 'YData', [0, Cy]);
    set(slider_D, 'XData', Dx, 'YData', Dy);
    set(slider_C, 'XData', Cx, 'YData', Cy);

    pause(0.01); % Control animation speed
end
