% Simulation of Differential Drive Robot with Two-Link Manipulator using ODE45
clc; clear; close all;

% Parameters
m = 10;         % Mass of the robot (kg)
I = 5;          % Moment of inertia of the robot (kg.m^2)
R = 0.5;        % Distance between wheels (m)
L1 = 1;         % Length of link 1 (m)
L2 = 0.8;       % Length of link 2 (m)
g = 9.81;       % Gravity (m/s^2)
dt = 0.1;       % Time step (s)
T = 10;         % Total simulation time (s)
time = 0:dt:T;  % Time vector

% Initial state vector: [x; y; theta; q1; q2; v; omega; dq1; dq2]
initial_state = [0; 0; 0; 0; 0; 0; 0; 0; 0];

% Solve ODE using ode45
[t, state_history] = ode45(@(t, state) system_dynamics(t, state, m, I, R, L1, L2, g), time, initial_state);

% Extract states for plotting and animation
x = state_history(:, 1);
y = state_history(:, 2);
theta = state_history(:, 3);
q1 = state_history(:, 4);
q2 = state_history(:, 5);

%% Animation
animate_robot(x, y, theta, q1, q2, L1, L2, t);

% System dynamics function
function dstate = system_dynamics(~, state, m, I, R, L1, L2, g)
    % Unpack state vector
    x = state(1);
    y = state(2);
    theta = state(3);
    q1 = state(4);
    q2 = state(5);
    v = state(6);
    omega = state(7);
    dq1 = state(8);
    dq2 = state(9);

    % Random control inputs
    Fr = 5*randn;         % Force on right wheel
    Fl = 8*randn;         % Force on left wheel
    tau1 = 6*randn;       % Torque on joint 1
    tau2 = 9*randn;       % Torque on joint 2

    % Dynamics of the differential drive robot
    dv = (Fr + Fl) / m;                % Linear acceleration
    domega = R * (Fr - Fl) / I;        % Angular acceleration

    % Dynamics of the manipulator (simplified)
    M = [1 0; 0 1];                    % Inertia matrix (identity for simplicity)
    C = [0 0; 0 0];                    % Coriolis matrix (neglected)
    G = [0; 0];                        % Gravity vector (neglected)
    tau = [tau1; tau2];                % Joint torques
    ddq = M \ (tau - C * [dq1; dq2] - G); % Joint accelerations

    % State derivatives
    dx = v * cos(theta);
    dy = v * sin(theta);
    dtheta = omega;
    dq1 = dq1;
    dq2 = dq2;
    dv = dv;
    domega = domega;
    ddq1 = ddq(1);
    ddq2 = ddq(2);

    % Pack state derivatives
    dstate = [dx; dy; dtheta; dq1; dq2; dv; domega; ddq1; ddq2];
end

% Animation function
function animate_robot(x, y, theta, q1, q2, L1, L2, t)
    figure;
    axis equal;
    grid on;
    xlabel('x (m)');
    ylabel('y (m)');
    title('Robot and Manipulator Animation');
    xlim([min(x) - 2, max(x) + 2]);
    ylim([min(y) - 2, max(y) + 2]);

    for i = 1:length(t)
        % Clear previous frame
        cla;

        % Robot position and orientation
        robot_x = x(i);
        robot_y = y(i);
        robot_theta = theta(i);

        % Draw robot body
        rectangle('Position', [robot_x - 0.2, robot_y - 0.1, 0.4, 0.2], 'Curvature', 0.2, 'FaceColor', [0.8, 0.8, 1]);
        hold on;

        % Draw robot wheels
        wheel1_x = robot_x - 0.2 * cos(robot_theta);
        wheel1_y = robot_y - 0.2 * sin(robot_theta);
        wheel2_x = robot_x + 0.2 * cos(robot_theta);
        wheel2_y = robot_y + 0.2 * sin(robot_theta);
        plot([wheel1_x, wheel2_x], [wheel1_y, wheel2_y], 'k', 'LineWidth', 2);

        % Draw manipulator
        joint1_x = robot_x + L1 * cos(q1(i) + robot_theta);
        joint1_y = robot_y + L1 * sin(q1(i) + robot_theta);
        joint2_x = joint1_x + L2 * cos(q1(i) + q2(i) + robot_theta);
        joint2_y = joint1_y + L2 * sin(q1(i) + q2(i) + robot_theta);

        plot([robot_x, joint1_x], [robot_y, joint1_y], 'r', 'LineWidth', 2);
        plot([joint1_x, joint2_x], [joint1_y, joint2_y], 'b', 'LineWidth', 2);

        % Pause for animation
        pause(0.05);
    end
end