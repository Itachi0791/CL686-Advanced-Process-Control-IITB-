clc; clear; close all;

%% System Parameters
m = 5;      % Mass of hopper (kg)
k = 1000;   % Spring constant (N/m)
c = 0;      % Damping coefficient (N.s/m)
g = 9.81;   % Gravity (m/s^2)
l0 = 1.0;   % Natural leg length (m)

Ts = 0.01;   % Sampling time
N = 10;      % Prediction horizon

%% Define State Space Model
% State: [x, x_dot, y, y_dot, l, l_dot, theta, theta_dot]
A = [  0  1  0  0  0  0  0  0;
        0  0  0  0  0 -k/m 0  0;
        0  0  0  1  0  0  0  0;
        0  0  0  0  0 -k/m 0  0;
        0  0  0  0  0  1  0  0;
        0  0  0  0 -k/m -c/m 0  0;
        0  0  0  0  0  0  0  1;
        0  0  0  0  0  0  0  0];
B = [0 0;
     1/m 0;
     0 0;
     0 1/m;
     0 0;
     0 0;
     0 0;
     0 0];

%% MPC Weights and Constraints
Q = diag([10 10 10 10 10 10 10 10]);
R = diag([0.1 0.1]);
umin = [-pi/4; 0];
umax = [pi/4; 200];

%% Simulation Setup
Tsim = 5; % Total simulation time (seconds)
time = 0:Ts:Tsim;
x = zeros(8, length(time)); % State trajectory
x(:,1) = [0; 0; 1; 0; 1; 0; 0; 0]; % Initial conditions
u = zeros(2, length(time));

% Desired outputs
xd_ref = 2; % Desired forward speed (m/s)
y_ref = 1.2; % Desired hopping height (m)

%% MPC Optimization Loop
for k = 1:length(time)-1
    % Define reference trajectory
    ref = [0; xd_ref; y_ref; 0; l0; 0; 0; 0];
    
    % Solve quadratic program for optimal control
    H = 2 * (R + B' * Q * B);
    f = 2 * (x(:,k)' * A' * Q * B - ref' * Q * B);
    u = quadprog(H, f, [], [], [], [], umin, umax);
    
    % Apply control and integrate dynamics
    x(:,k+1) = A*x(:,k) + B*u;
    
    % Store control inputs
    u(:,k) = u;
end

%% Animation Setup
figure;
for k = 1:10:length(time)
    clf;
    hold on; axis equal;
    x_hop = x(1,k);
    y_hop = x(3,k);
    l_hop = x(5,k);
    theta_hop = x(7,k);
    
    % Draw ground
    plot([-5 5], [0 0], 'k', 'LineWidth', 2);
    
    % Draw leg
    leg_x = [x_hop x_hop - l_hop * sin(theta_hop)];
    leg_y = [y_hop y_hop - l_hop * cos(theta_hop)];
    plot(leg_x, leg_y, 'b', 'LineWidth', 2);
    
    % Draw body
    viscircles([x_hop, y_hop], 0.1, 'Color', 'r');
    
    pause(0.01);
end

%% Plot Results
figure;
subplot(3,1,1);
plot(time, x(2,:)); hold on;
yline(xd_ref, 'r--');
title('Forward Speed Regulation'); xlabel('Time (s)'); ylabel('x-dot (m/s)');

subplot(3,1,2);
plot(time, x(3,:)); hold on;
yline(y_ref, 'r--');
title('Hopping Height Regulation'); xlabel('Time (s)'); ylabel('y (m)');

subplot(3,1,3);
plot(time(1:end-1), u(1,:)); hold on;
plot(time(1:end-1), u(2,:));
title('Control Inputs'); xlabel('Time (s)'); ylabel('Inputs'); legend('Theta','Thrust');
