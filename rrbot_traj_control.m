clear; close; clc;
syms m1 m2 theta1 theta2 r1 r2 l1 l2 I1 I2 theta1_dot theta1_ddot theta2_dot theta2_ddot u1 u2 g 'real'
m1=1; m2=1; l1=1; l2=1 ;r1=0.45; r2=0.45; g=9.81 ;I2= 0.084; I1= 0.084;

% ROS Setup
rosinit;

j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');

JointStates = rossubscriber('/rrbot/joint_states');

tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);

tic;
t = 0;
i = 1;

while(t < 10)
    t = toc;

    % read the joint states
    jointData = receive(JointStates);
    x = [jointData.Position(1); jointData.Position(2);
      jointData.Velocity(1); jointData.Velocity(2)];

    q1_d = pi - (3*pi.*t.^2)/100 + (pi.*t.^3)/500;
    q1dot_d = - (6*pi.*t)/100 + (3.*t.^2*pi)/500;
    q1ddot_d = - (6*pi)/100 + (6.*t*pi)/500;
    
    q2_d =  pi/2 - (3*pi.*t.^2)/200 + (pi.*t.^3)/1000;
    q2dot_d = - (6*pi.*t)/200 + (3.*t.^2*pi)/1000;
    q2ddot_d = - (6*pi)/200 + (6.*t*pi)/1000;

    theta1 = x(1);
    theta2 = x(2);
    theta1_dot = x(3);
    theta2_dot = x(4);

    M = [m2*l1^2 + 2*m2*cos(theta2)*l1*r2 + m1*r1^2 + m2*r2^2 + I1 + I2, m2*r2^2 + l1*m2*cos(theta2)*r2 + I2; I2 + m2*r2^2+ m2*r2*l1*cos(theta2), m2*r2^2 + I2] ;
    C1 = [-theta2_dot*((9*sin(theta2)*(theta1_dot + theta2_dot))/20 + (9*theta1_dot*sin(theta2))/20)];
    C2 = [(9*theta1_dot*sin(theta2)*(theta1_dot + theta2_dot))/20 - (9*theta1_dot*theta2_dot*sin(theta2))/20];
    G = [- g*l1*m2*sin(theta1) - g*m1*r1*sin(theta1) - m2*g*r2*sin(theta1 + theta2); - g*m2*r2*sin(theta1 + theta2)];

    k = [110, 0 , 21, 0; 0, 110, 0, 21]; 
    e = x - [q1_d; q2_d; q1dot_d; q2dot_d];
    e(1:2, 1) = wrapToPi(e(1:2,1));   

    % State feedback controller
    U = M*(- k*(e) + [q1ddot_d; q2ddot_d]) + [C1; C2] + G;

    tau1.Data = [U(1,:)];
    tau2.Data = [U(2,:)];
    
    send(j1_effort, tau1);
    send(j2_effort, tau2);

    % sampling data to be plotted
    g1(i) = x(1);
    g2(i) = x(2);
    g3(i) = x(3);
    g4(i) = x(4);

    q1d(i) = pi - (3*pi.*t.^2)/100 + (pi.*t.^3)/500;
    q1dotd(i) = - (6*pi.*t)/100 + (3.*t.^2*pi)/500;
    
    q2d(i) =  pi/2 - (3*pi.*t.^2)/200 + (pi.*t.^3)/1000;
    q2dotd(i) = - (6*pi.*t)/200 + (3.*t.^2*pi)/1000;
    
    force1(i) = tau1.Data;
    force2(i) = tau2.Data;
    
    time(i) = t;
    
    i = i + 1;
end

figure(1)
subplot(2, 2, 1)
plot(time, q1d);
xlabel('t', 'FontSize',14)
ylabel('q1_d','FontSize',14);
hold on
plot(time,g1,'r');
xlabel('t', 'FontSize',14)
ylabel('theta1','FontSize',14);
legend('desired','actual')
hold off

subplot(2, 2, 2)
plot(time, q1dotd);
xlabel('t', 'FontSize',14)
ylabel('q1_d','FontSize',14);
hold on
plot(time,g3,'r');
xlabel('t', 'FontSize',14)
ylabel('theta1 dot','FontSize',14)
legend('desired','actual')
hold off

subplot(2, 2, 3)
plot(time, q2d);
xlabel('t', 'FontSize',14)
ylabel('q2_d','FontSize',14);
hold on
plot(time,g2,'r');
xlabel('t', 'FontSize',14)
ylabel('theta2','FontSize',14)
legend('desired','actual')
hold off

subplot(2, 2, 4)
plot(time, q2dotd);
xlabel('t', 'FontSize',14)
ylabel('q2_d','FontSize',14);
hold on
plot(time,g4,'r');
xlabel('t', 'FontSize',14)
ylabel('theta2 dot','FontSize',14)
legend('desired','actual')
hold off


figure(2)
subplot(2,1,1);
plot(time,force1);
xlabel('t');
ylabel('u1');

figure(2)
subplot(2,1,2);
plot(time,force2);
xlabel('t');
ylabel('u2');

tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);

% disconnect from roscore
rosshutdown;