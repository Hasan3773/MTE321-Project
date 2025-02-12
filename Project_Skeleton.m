%% Project Skeleton Code

clear; clc; close all;

%%initial parameter: unit: m, degree, rad/sec
r2 = 36; % cm  o2B
% and so on ...


theta2 = 0:1:360; % from 0 to 360 with step 1: [0,1,2,3,4....360]
dtheta2 = 1;
ddtheta2 = 0; 

% TIPS:  

% cosd(x) - is a cosine of x, where x in degrees
% cos(x) - is a cosine of x, where x in radians
% using '.*' enables element-wise multiplication
% accordingly, '.^' element-wise exponent
% [a1 a2 a3].^[b1 b2 b3] = [a1*b1 a2*b2 a3*b3]
% '*' is matrix multiplication

%% Part 1- Calculations for kinematic variables, caculated based on loop closure eqn

r_i = % ENTER YOUR CODE HERE %
theta_i = % ENTER YOUR CODE HERE %
% Hint: Check if the angle needs to be adjusted to its true value
% Hint: Check this for all other angles too



%% Take time derivative of loop eqn (d/dt) 
% and solve them for dtheta3, dtheta5 & dr6
% and the same for the second derivatives. 

dr_i = % ENTER YOUR CODE HERE %
dtheta_i = % ENTER YOUR CODE HERE %

ddr_i = % ENTER YOUR CODE HERE %
ddtheta_i = % ENTER YOUR CODE HERE %

% and so on


%% Plot vars;

% Plot all desired deliverables. 

figure (1)
plot(theta_2,theta_i)
grid on;
title('$\theta_i$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('\theta_i   unit: degree')

% and so on

% *****************************************************

%% Part 2 - Force and Moment Calculation


%%initial parameters:

dtheta_2 = 2; % theta2 dot
ddtheta_2 = 0; % theta2 doble-dot - second derivative

rho = % ENTER YOUR CODE HERE %; % density, gr/cm3
d = % ENTER YOUR CODE HERE %; % diameter, cm

m_i = % ENTER YOUR CODE HERE % ; % link 2, o2a2 kg
I_Gi = % ENTER YOUR CODE HERE %;
% and so on


M12_list = [];
theta2_list = [];
Fs_list = [];  % shaking force
alpha_s_list = []; % direction of a shaking force
Ms_list =[]; % Shaking moment
Fij_list = []; % Forces
Fij_alpha = []; % Angles at which forces are acting


for theta2 = 0:1:360

    % kinematic variables are caculated based on loop eqn
    r_i = % ENTER YOUR CODE HERE %;
    dr_i = % ENTER YOUR CODE HERE %;
    ddr_i = % ENTER YOUR CODE HERE %;

% and so on    

    B = get_ma_vector(%m_i, ... % these are the examples of the possible input
        % ri ... % Only include the inputs that are necessary
        % theta_i ...
        % dtheta_i ...
        % ddtheta_i ...
        % ddr_i, ...
        % I_Gi);
    
    A = get_A_matrix(%m_i, ... % these are the examples of the possible input
        % ri ... % Only include the inputs that are necessary
        % theta_i ...
        % dtheta_i ...
        % ddtheta_i ...
        % ddr_i, ...
        % I_Gi);

    x = A\ B; % Ax = B, solution for x; note that in MATLAB: A\B = B/A
    
    % M12:
    M12 = x(% ENTER YOUR CODE HERE%);
    M12_list = [M12_list; M12];
    
    Fijx = x(% ENTER YOUR CODE HERE%);
    Fijy = x(% ENTER YOUR CODE HERE%);
    
    % Magnitudes of all forces: 
    % Atan is defined on [-pi/2; pi/2]. 
    % This if clause will help to adjust the value of the angle 
    % to its true value:	
    Fij_list = [Fij_list; % ENTER YOUR CODE HERE%];

    
    % Directions of all forces:    
    fx = % ENTER YOUR CODE HERE%;
    fy = % ENTER YOUR CODE HERE%;
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end 
    Fij_alpha = [Fij_alpha; alpha_f];

    % and so on
    
  
    % Collecting the values of theta2:
    theta2_list = [theta2_list, theta_2];
     
   
    
end


% Regular and Polar plots:
% Might have to transpose the Force vectors for polar plot. Do so if needed
% Polar plot only works with radians so will have to do it accordingly

figure (3)
plot(theta2_list,M12_list)
grid on;
title('M_{12} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('M12   unit: N-m')


% Convert degrees to the radians
theta2_rad = deg2rad(theta2_list);

%figure (4)
%polarplot(Fij_alpha,Fij_list)
%grid on;
%title('F_{ij} polar plot')

% and so on ...