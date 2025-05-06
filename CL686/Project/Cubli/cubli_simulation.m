clear; clc; close all;

% Define system parameters
params.mh = 0.5;  % Mass of Cubli (kg)
params.mw = 0.1;  % Mass of wheels (kg)
params.l = 0.05;  % Distance from pivot to CoM (m)
params.g = 9.81;  % Gravity (m/s^2)
params.Jh = 0.02; % Moment of inertia of Cubli (kg·m²)
params.Jw = 0.005; % Moment of inertia of wheels (kg·m²)
params.Km = 0.05; % Motor constant
params.Cw = 0.01; % Damping constant

% Initial conditions (small perturbation)
theta0 = 5 * pi/180; % 5-degree tilt
omega0 = 0;
omega_w0 = 0;
state0 = [theta0; omega0; omega_w0]; % [angle; angular velocity; wheel velocity]

% Time span
tspan = [0 10];

% Simulate Nonlinear System
[t, y] = ode45(@(t, x) cubli_nonlinear(t, x, params), tspan, state0);

% Simulate with LQR Control
[K, A, B] = cubli_lqr(params);
[t2, y2] = ode45(@(t, x) cubli_controlled(t, x, K, params), tspan, state0);

% Plot results
figure;
subplot(2,1,1);
plot(t, y(:,1) * 180/pi, 'r', 'LineWidth', 1.5); hold on;
plot(t2, y2(:,1) * 180/pi, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angle (degrees)');
title('Cubli Balance: Nonlinear vs LQR Control');
legend('Uncontrolled', 'LQR Controlled');
grid on;

subplot(2,1,2);
plot(t, y(:,2), 'r', 'LineWidth', 1.5); hold on;
plot(t2, y2(:,2), 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocity');
grid on;
