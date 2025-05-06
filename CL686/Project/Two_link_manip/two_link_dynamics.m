clc; clear; close all;

% Simulation parameters
dt = 0.02; t_end = 10; tspan = 0:dt:t_end;

% Physical parameters
m1 = 5; m2 = 4; l1 = 1; l2 = 0.8; I1 = m1*l1^2/3; I2 = m2*l2^2/3; g = 0;
u = [3; -3]; % Control input (torque)

% Initial conditions [theta1, theta2, theta1_dot, theta2_dot]
q0 = [pi/2; pi/4; 0; 0];

% Dynamics matrices
H = @(q) [I1 + m1*l1^2/4 + I2 + m2*(l2^2/4 + l1^2 + l1*l2*cos(q(2))), I2 + m2*(l2^2/4 + l1*l2*cos(q(2))/2); 
          I2 + m2*(l2^2/4 + l1*l2*cos(q(2))/2), I2 + m2*l2^2/4];
C = @(q) [-m2*l1*l2*sin(q(2))*q(4), -m2*l1*l2*sin(q(2))*q(4)/2; 
           m2*l1*l2*sin(q(2))*q(3)/2, 0];
G = @(q) [m1*g*l1*cos(q(1))/2 + m2*g*(l1*cos(q(1)) + l2/2*cos(q(1)+q(2))); 
          m2*g*l2/2*cos(q(1)+q(2))];

% Dynamics function for ode89
arm_dynamics = @(t, q) [q(3:4); H(q) \ (u - C(q)*q(3:4) - G(q))];

% Integrate using ode89
options = odeset('abstol', 1e-9, 'reltol', 1e-9);
[t_out, q_out] = ode89(arm_dynamics, tspan, q0, options);

% Compute end-effector trajectory
trajectory = [l1*cos(q_out(:,1)) + l2*cos(q_out(:,1)+q_out(:,2)), ...
              l1*sin(q_out(:,1)) + l2*sin(q_out(:,1)+q_out(:,2))];

% Animation
figure; axis equal; grid on; hold on;
xlim([-1.1*(l1+l2), 1.1*(l1+l2)]); ylim([-1.1*(l1+l2), 1.1*(l1+l2)]);
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('Manipulator Animation with End-Effector Trajectory');

% Initialize plot handles
h_links = plot([0, 0], [0, 0], 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 2);
h_trajectory = plot(trajectory(:,1), trajectory(:,2), 'b-', 'LineWidth', 1.5);
h_end_effector = plot(0, 0, 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2);

for i = 1:length(t_out)
    q = q_out(i, :);
    x1 = l1*cos(q(1)); y1 = l1*sin(q(1));
    x2 = x1 + l2*cos(q(1)+q(2)); y2 = y1 + l2*sin(q(1)+q(2));
    
    % Update manipulator links
    set(h_links, 'XData', [0, x1, x2], 'YData', [0, y1, y2]);
    
    % Update end-effector position
    set(h_end_effector, 'XData', x2, 'YData', y2);
    
    % Update trajectory (up to current point)
    set(h_trajectory, 'XData', trajectory(1:i,1), 'YData', trajectory(1:i,2));
    
    pause(dt);
end