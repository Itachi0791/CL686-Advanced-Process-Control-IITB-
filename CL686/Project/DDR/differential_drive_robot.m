clc; clear; close all;

%% Simulation Parameters
dt = 0.1;  % Time step (s)
T = 20;    % Total simulation time (s)
tspan = 0:dt:T;

%% Robot Parameters
wheel_radius = 0.05;  % Wheel radius (m)
wheel_base = 0.3;     % Distance between wheels (m)

% Initial State [x, y, theta]
state0 = [0; 0; 0]; % (x, y, heading)

% Velocity Commands (v, omega)
v = 2;  % Linear velocity (m/s)
w = 1;  % Angular velocity (rad/s)

%% Solve Using ode45
[t, trajectory] = ode45(@(t, state) diff_drive_ode(t, state, v, w), tspan, state0);

%% Plot Simulation
figure; hold on; axis equal;
xlabel('X (m)'); ylabel('Y (m)');
title('Differential Drive Robot Simulation');

for i = 1:length(t)
    clf; hold on; axis equal;
    plot(trajectory(1:i, 1), trajectory(1:i, 2), 'b', 'LineWidth', 2);
    plot_robot(trajectory(i, 1), trajectory(i, 2), trajectory(i, 3), wheel_base);
    xlabel('X (m)'); ylabel('Y (m)');
    title('Differential Drive Robot Simulation');
    pause(0.05);
end

%% Function for Differential Drive Model
function dstate = diff_drive_ode(~, state, v, w)
    x = state(1);
    y = state(2);
    theta = state(3);
    
    dx = v * cos(theta);
    dy = v * sin(theta);
    dtheta = w;
    
    dstate = [dx; dy; dtheta];
end

%% Function to Draw Robot
function plot_robot(x, y, theta, wheel_base)
    % Robot body
    L = 0.4; W = 0.2;
    body = [-L/2, L/2, L/2, -L/2, -L/2;
            -W/2, -W/2, W/2, W/2, -W/2];
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    body = R * body;
    fill(body(1,:) + x, body(2,:) + y, 'r');

    % Wheels
    wheel1 = [-wheel_base/2, wheel_base/2;
              -W/2, -W/2];
    wheel2 = [-wheel_base/2, wheel_base/2;
               W/2,  W/2];
    wheel1 = R * wheel1;
    wheel2 = R * wheel2;
    plot(wheel1(1,:) + x, wheel1(2,:) + y, 'k', 'LineWidth', 3);
    plot(wheel2(1,:) + x, wheel2(2,:) + y, 'k', 'LineWidth', 3);
end

