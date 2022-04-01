clear; clc; close all;
syms m1 m2 theta1 theta2 r1 r2 l1 l2 I1 I2 theta1_dot theta1_ddot theta2_dot theta2_ddot u1 u2 g 'real'
syms a0 a1 a2 a3 t to tf qo qodot qf qfdot
m1=1; m2=1; l1=1; l2=1 ;r1=0.45; r2=0.45; g=9.81 ;I2= 0.084; I1= 0.084;

% Generation of a cubic trajectory
A = [1 to to^2 to^3; 0 1 2*to 3*(to^2); 1 tf tf^2 tf^3; 0 1 2*tf 3*(tf^2)];
b = [qo; qodot; qf; qfdot];
a = inv(A) * b;

% Coefficient matrix of the first joint
a1 = subs(a, [to, tf, qo, qodot, qf, qfdot], [0, 10, pi, 0, 0 0]);

% Coefficient matrix of the second joint
a2 = subs(a, [to, tf, qo, qodot, qf, qfdot], [0, 10, pi/2, 0, 0 0]);

cubic_trajec1 = a1(1) + a1(2)*t + a1(3)*(t^2) + a1(4)*(t^3);
cubic_trajec2 = a2(1) + a2(2)*t + a2(3)*(t^2) + a2(4)*(t^3);

eq1 = theta1_ddot*(m2*l1^2 + 2*m2*cos(theta2)*l1*r2 + m1*r1^2 + m2*r2^2 + I1 + I2) - theta2_dot*(l1*m2*r2*theta1_dot*sin(theta2) + l1*m2*r2*sin(theta2)*(theta1_dot + theta2_dot)) - u1 + theta2_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) - g*l1*m1*sin(theta1) - g*m1*r1*sin(theta1);
eq2 = theta2_ddot*(m2*r2^2 + I2) - u2 + theta1_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) + l1*m2*r2*theta1_dot*sin(theta2)*(theta1_dot + theta2_dot) - l1*m2*r2*theta1_dot*theta2_dot*sin(theta2);

G1 = subs(eq1, [theta1_ddot, theta1_dot, theta2_ddot, theta2_dot], [0, 0, 0, 0]);
G2 = subs(eq2, [theta1_ddot, theta1_dot, theta2_ddot, theta2_dot], [0, 0, 0, 0]);

M1 = subs(eq1, [theta1_dot, theta2_dot], [0, 0]);
M2 = subs(eq2, [theta1_dot, theta2_dot], [0, 0]);

C1 = subs(eq1, [theta1_ddot, theta2_ddot], [0, 0]) - G1;
C2 = subs(eq2, [theta1_ddot, theta2_ddot], [0, 0]) - G2;

M = [m2*l1^2 + 2*m2*cos(theta2)*l2*r2 + m1*r1^2 + m2*r2^2 + I1 + I2, m2*r2^2 + l1*m2*cos(theta2)*r2 + I2; I2 + m2*r2^2+ m2*r2*l1*cos(theta2), m2*r2^2 + I2] ;
C= [l1*m2*r2*sin(theta2), -l1*m2*r2*theta1_dot*sin(theta2) + l1*m2*r2*sin(theta2); l1*m2*r2*sin(theta2)*(theta1_dot + theta2_dot) - l1*m2*r2*theta2_dot*sin(theta2), 0];
G = [- g*l1*m2*sin(theta1) - g*m1*r1*sin(theta1) - m2*g*r2*sin(theta1+theta2); - g*m2*r2*sin(theta1 + theta2)];

% simplify(M*[theta1_ddot; theta2_ddot] + C*[theta1_dot; theta2_dot] + G - [eq1; eq2])
% simplify([C1; C2] - C*[theta1_dot; theta2_dot])

% Eigenvalue placement method
lambda = [-1, -2, -3, -4];
A = [0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 0, 0; 0, 0, 0, 0];
B = [0, 0; 0, 0; 1, 0; 0, 1];
K = place(A, B, lambda)

% Simulation of the system 
T = 10;
[t,y] = ode45(@ode_link, [0,T], [deg2rad(200),deg2rad(125),0, 0]);
% [t,y] = ode45(@ode_link, [0,T], [pi, pi/2, 0, 0]); % Start from the desired trajectory's initial conditions
% [t,y] = ode45(@ode_link, [0,T], [0, 0, 0, 0]); % To start from 0

K = [12, 0 , 7, 0; 0, 2, 0, 3];

q1_d = pi - (3*pi.*t.^2)/100 + (pi.*t.^3)/500;
q1dot_d = - (6*pi.*t)/100 + (3.*t.^2*pi)/500;
q1ddot_d = - (6*pi)/100 + (6.*t*pi)/500;

q2_d =  pi/2 - (3*pi.*t.^2)/200 + (pi.*t.^3)/1000;
q2dot_d = - (6*pi.*t)/200 + (3.*t.^2*pi)/1000;
q2ddot_d =- (6*pi)/200 + (6.*t*pi)/1000;

for i = 1:size(y,1)
    v1(i) = - K(1,1)*(y(i,1) - q1_d(i)) - K(1,2)*(y(i,2) - q2_d(i)) - K(1,3)*(y(i,3) - q1dot_d(i)) - K(1,4)*(y(i,4) - q1dot_d(i)) + q1ddot_d(i);
    v2(i) = - K(2,1)*(y(i,1) - q1_d(i)) - K(2,2)*(y(i,2) - q2_d(i)) - K(2,3)*(y(i,3) - q1dot_d(i)) - K(2,4)*(y(i,4) - q1dot_d(i)) + q2ddot_d(i);
    
    M = [m2*l1^2 + 2*m2*cos(y(i,2))*l2*r2 + m1*r1^2 + m2*r2^2 + I1 + I2, m2*r2^2 + l1*m2*cos(y(i,2))*r2 + I2; I2 + m2*r2^2+ m2*r2*l1*cos(y(i,2)), m2*r2^2 + I2] ;
    C = [-y(i,4)*((9*sin(y(i,2))*(y(i,3) + y(i,4)))/20 + (9*y(i,3)*sin(y(i,2)))/20); (9*y(i,3)*sin(y(i,2))*(y(i,3) + y(i,4)))/20 - (9*y(i,3)*y(i,4)*sin(y(i,2)))/20];
    G = [- g*l1*m2*sin(y(i,1)) - g*m1*r1*sin(y(i,1)) - m2*g*r2*sin(y(i,1) + y(i,2)); - g*m2*r2*sin(y(i,1) + y(i,2))];

    u1(i)= M(1, :) * [v1(i); v2(i)] + C(1) + G(1) ;
    u2(i)= M(2, :) * [v1(i); v2(i)] + C(2) + G(2) ;
end

figure(1)
subplot(2,2,1);
plot(t, q1_d);
xlabel('t', 'FontSize',14)
ylabel('q1_d','FontSize',14);
hold on
plot(t,y(:,1),'r');
xlabel('t', 'FontSize',14)
ylabel('theta1','FontSize',14);
legend('desired','actual')
hold off

subplot(2,2,2);
plot(t, q1dot_d);
xlabel('t', 'FontSize',14)
ylabel('q1_d','FontSize',14);
hold on
plot(t,y(:,3),'r');
xlabel('t', 'FontSize',14)
ylabel('theta1 dot','FontSize',14)
legend('desired','actual')
hold off

subplot(2,2,3);
plot(t, q2_d);
xlabel('t', 'FontSize',14)
ylabel('q2_d','FontSize',14);
hold on
plot(t,y(:,2),'r');
xlabel('t', 'FontSize',14)
ylabel('theta2','FontSize',14)
legend('desired','actual')
hold off

subplot(2,2,4);
plot(t, q2dot_d);
xlabel('t', 'FontSize',14)
ylabel('q2_d','FontSize',14);
hold on
plot(t,y(:,4),'r');
xlabel('t', 'FontSize',14)
ylabel('theta2 dot','FontSize',14)
legend('desired','actual')
hold off

figure(2)
subplot(2, 1, 1);
plot(t, u1,'b');
xlabel('t', 'FontSize', 14)
ylabel('u1', 'FontSize', 14);

subplot(2, 1, 2);
plot(t, u2,'b');
xlabel('t', 'FontSize', 14)
ylabel('u2', 'FontSize', 14);

function dX = ode_link(t,X)
m1=1;m2=1; l1=1; l2=1 ;r1=0.45; r2=0.45; g=9.81 ;I2= 0.084; I1= 0.084;

dX= zeros(4,1);
X=num2cell(X);
[theta1, theta2, theta1_dot, theta2_dot] = deal(X{:});

M = [m2*l1^2 + 2*m2*cos(theta2)*l1*r2 + m1*r1^2 + m2*r2^2 + I1 + I2, m2*r2^2 + l1*m2*cos(theta2)*r2 + I2; I2 + m2*r2^2+ m2*r2*l1*cos(theta2), m2*r2^2 + I2] ;
C = [-theta2_dot*((9*sin(theta2)*(theta1_dot + theta2_dot))/20 + (9*theta1_dot*sin(theta2))/20); (9*theta1_dot*sin(theta2)*(theta1_dot + theta2_dot))/20 - (9*theta1_dot*theta2_dot*sin(theta2))/20];
G = [- g*l1*m2*sin(theta1) - g*m1*r1*sin(theta1) - m2*g*r2*sin(theta1 + theta2); - g*m2*r2*sin(theta1 + theta2)];

q1_d = pi - (3*pi*t^2)/100 + (pi*t^3)/500;
q1dot_d = - (6*pi*t)/100 + (3*t^2*pi)/500;
q1ddot_d = - (6*pi)/100 + (6*t*pi)/500;

q2_d =  pi/2 - (3*pi*t^2)/200 + (pi*t^3)/1000;
q2dot_d = - (6*pi*t)/200 + (3*t^2*pi)/1000;
q2ddot_d = - (6*pi)/200 + (6*t*pi)/1000;

% Feedback linearization control 
K = [12, 0 , 7, 0; 0, 2, 0, 3];

U = M*(- K*([theta1; theta2; theta1_dot; theta2_dot] - [q1_d; q2_d; q1dot_d; q2dot_d]) +[q1ddot_d; q2ddot_d]) + C + G;

u1 = [U(1,:)];
u2 = [U(2,:)];

dX(1) = theta1_dot;

dX(2) = theta2_dot;

dX(3) = (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

% ddq = M\(U - C * [theta1_dot; theta2_dot] - G);
% dX(3) = ddq(1, 1);
% dX(4) = ddq(2, 1);

dX(4) = -(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + g*l1^2*m1*m2*r2*cos(theta2)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

end

