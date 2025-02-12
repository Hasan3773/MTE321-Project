%% Animate the motion of the mechanism 
%%initial parameter: unit: cm, degree, rad/sec
r2 = 36; % cm  AB
% and so on

close all
inp_theta2 = 0:1:720; % two rotations from 0 to 360 with step 1: [0,1,2,3,4....360]

for i = 1:length(inp_theta2)

    % Loop-1

    % Define angles and displacements for loop1 
    theta_2 = i; 
    theta_i = % ENTER YOUR CODE HERE
    % Hint: Check if the angle needs to be adjusted to its true value
    % Hint: Check this for all other angles too
    r_i = % ENTER YOUR CODE HERE

   
    % Loop-1 (o2 is origin)
    % Plot Links
    % ENTER YOUR CODE HERE %;
    link_i = % ENTER YOUR CODE HERE %;
    
    % and so on    



    % Loop-2 
    
    % ENTER YOUR CODE HERE %;

 
    % Interval to update plot 
    pause(0.001);
    delete(link_i); 
    
    % and so on

    grid on;
   
end 
