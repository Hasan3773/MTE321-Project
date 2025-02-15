%% Animate the motion of the mechanism 
%%initial parameter: unit: cm, degree, rad/sec
r2 = 36; % cm  AB
r3 = 120; % variable
r6 = 60;
r8 = 8.4; 

close all
inp_theta2 = 0:1:720; % two rotations from 0 to 360 with step 1: [0,1,2,3,4....360]

for i = 1:length(inp_theta2)

    % Loop-1
    % Hint: Check if the angle needs to be adjusted to its true value
    % Hint: Check this for all other angles too
    % Define angles and displacements for loop1 
    
    % plot fixed pivots
    figure(1)
    plot(0,0,'o')
    hold on

    theta_2 = inp_theta2(i); 
    Ax = 0; Ay = 0; % fixed pivot
    Bx = r2 * cosd(theta_2);
    By = r2 * sind(theta_2);

    theta_3 = asind((r2*sind(theta_2) - r8) ./ r3);
    theta_3 = theta_3 + 180; 
    theta_5 = theta_3;
    r7 = (r2 * cosd(theta_2)) - (r3 * cosd(theta_3)); 
    theta_6 = 180 - asind((r2*sind(theta_2 + theta_5) / r6)) - theta_5;
    
    r5 = ((-r2 * sind(theta_2 + theta_6)) / sind(theta_5 + theta_6));

    Dy = r8; % constraint
    Dx = r7;
    plot(Dx,Dy,'o')

    Cx = r6 * cosd(theta_6);
    Cy = r6 * sind(theta_6);

    % crank (A -> B)
    link_AB = plot([Ax, Bx], [Ay, By], 'ro-', 'LineWidth', 2);
    
    % Coupler (B -> D)
    link_BD = plot([Bx, Dx], [By, Dy], 'bo-', 'LineWidth', 2);
    
    % Ground link (A -> C)
    link_AC = plot([Ax, Cx], [Ay, Cy], 'ko-', 'LineWidth', 2);

    % Ground link (B -> C)
    % link_BC = plot([Ax, Cx], [Ay, Cy], 'ko-', 'LineWidth', 2);
    
    axis([-50 250 -50 250])
    pause(0.001);
    hold off;
    grid on;
   
end 
