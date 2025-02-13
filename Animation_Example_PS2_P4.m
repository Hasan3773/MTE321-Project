clear; clc; close all;

% define link lengths
r1=3.0;
r2=6.0;
r4=10.0;

% input angle theta2 from 0-360 deg
x = 0:1:360; 

% compute output angle theta3
theta3=atand((r1+r2.*sind(x))./(r2.*cosd(x)-r4)); % inverse tangent 
theta3=theta3+180; % make sure its in the right quadrant?

% calculate length of r3 ig
r3=sqrt(r1^2+r2^2+r4^2-2*r2*r4.*cosd(x)+2*r1*r2.*sind(x));

% plot output angle vs input angle
figure (1)
plot(x,theta3)
grid on;
title('$\theta_3$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('\theta_3   unit: degree')

% plot coupler length vs input angle 
figure (2)
plot(x,r3)
grid on;
title('$r_3$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('r_3   unit: cm')

for x = 0:1:720
theta3=atand((r1+r2.*sind(x))./(r2.*cosd(x)-r4));
theta3=theta3+180;
r3=sqrt(r1^2+r2^2+r4^2-2*r2*r4.*cosd(x)+2*r1*r2.*sind(x));

% plot fixed pivots
figure(3)
plot(0,0,'o')
hold on
plot(10,-3,'o')
hold on

% array to store [fixed pivot, crank endpoint, coupler endpoint]
A = [0 r2.*cosd(x) r2.*cosd(x)-18.*cosd(theta3)]; % x coord
B = [0 r2.*sind(x) r2.*sind(x)-18.*sind(theta3)]; % y coord
plot(0,0)
plot(A,B) % plot linkages
axis equal
axis([-6 26 -16 6])
line(A,B)
pause(0.001);
hold off

end

