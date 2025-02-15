function A = get_A_matrix(theta2,theta3,theta6)

% Matrix defined by equations underlined in blue in the assignment in order
%   I: 1    2    3                  4                 5                  6                  7    8    9                  10                11  12                13
% Var: F12x F12y F23x               F23y              F34x               F34y               F16x F16y F56x               F56y              F14 F35               M12
A = [  -1   0    -1                 0                 0                  0                  0    0    0                  0                 0   0                 0;
       0    -1   0                 -1                 0                  0                  0    0    0                  0                 0   0                 0;
       0    0    0.36*sind(theta2) -0.36*cosd(theta2) 0                  0                  0    0    0                  0                 0   0                 -1;
       0    0    -1                0                  -1                 0                  0    0    0                  0                 0   cosd(theta3+90)   0;
       0    0    0                 -1                 0                  -1                 0    0    0                  0                 0   sind(theta3+90)   0;
       0    0    0.60*sind(theta3) 0.60*cosd(theta3)  -0.60*sind(theta3) -0.60*cosd(theta3) 0    0    0                  0                 0   0                 0;
       0    0    0                 0                  -1                 0                  0    0    0                  0                 0   0                 0;
       0    0    0                 0                  0                  -1                 0    0    0                  0                 1   0                 0;
       0    0    0                 0                  0                  0                  0    0    1                  0                 0   -cosd(theta3-90)  0;
       0    0    0                 0                  0                  0                  0    0    0                  1                 0   -sind(theta3-90)  0;
       0    0    0                 0                  0                  0                  -1   0    1                  0                 0   0                 0;
       0    0    0                 0                  0                  0                  0    -1   0                  1                 0   0                 0;
       0    0    0                 0                  0                  0                  0    0    -0.60*sind(theta6) 0.60*cosd(theta6) 0   0                 0;
    ];


end