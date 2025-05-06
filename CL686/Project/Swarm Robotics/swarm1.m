clc; clear; close all;

%% Simulation Parameters
dt = 0.2;            % Time step
T = 10;              % Prediction horizon (number of steps)
N = 3;               % Number of robots
sim_time = 100;       % Simulation steps
v_max = 1.0;         % Max linear velocity
w_max = pi/4;        % Max angular velocity
d_min = 0.5;         % Minimum safe distance between robots

% Initial positions [x, y, theta]
robots = [0, 0, 0;
          2, 2, pi/4;
          -2, 1, -pi/4];

% Formation offsets (relative positions)
formation_offsets = [0, 0;    % Leader
                     1, 0.5;  % Right
                    -1, 0.5]; % Left

% Store trajectories
trajectories = cell(N,1);
for i = 1:N
    trajectories{i} = robots(i,1:2);
end

%% MPC Optimization Loop
for t = 1:sim_time
    % Define moving formation trajectory (leader moves along a sine wave)
    leader_x = 5 + 0.1 * t;          % Leader moves forward
    leader_y = 5 + sin(0.1 * t);     % Leader follows sine wave
    
    % Compute dynamic formation positions
    formation = [leader_x, leader_y;
                 leader_x + formation_offsets(2,1), leader_y + formation_offsets(2,2);
                 leader_x + formation_offsets(3,1), leader_y + formation_offsets(3,2)];
    
    for i = 1:N
        % Current state
        x0 = robots(i, 1);
        y0 = robots(i, 2);
        theta0 = robots(i, 3);
        
        % Set optimization variables: [v_1, w_1, v_2, w_2, ..., v_T, w_T]
        u0 = zeros(2*T, 1);  
        
        % Cost function
        cost_fun = @(u) cost_function(u, x0, y0, theta0, formation(i,:), T, dt);
        
        % Constraints
        A = []; b = [];
        Aeq = []; beq = [];
        lb = repmat([-v_max; -w_max], T, 1);
        ub = repmat([v_max; w_max], T, 1);
        
        % Solve optimization using fmincon
        options = optimoptions('fmincon', 'Display', 'none');
        u_opt = fmincon(cost_fun, u0, A, b, Aeq, beq, lb, ub, ...
                        @(u) collision_constraints(u, robots, i, T, dt, d_min), options);
        
        % Apply first control input
        v = u_opt(1);
        w = u_opt(2);
        robots(i, 1) = x0 + v * cos(theta0) * dt;
        robots(i, 2) = y0 + v * sin(theta0) * dt;
        robots(i, 3) = theta0 + w * dt;
        
        % Store trajectory
        trajectories{i} = [trajectories{i}; robots(i,1:2)];
    end

    % Plotting
    clf; hold on; grid on;
    for i = 1:N
        plot(trajectories{i}(:,1), trajectories{i}(:,2), 'b');
        plot(robots(i,1), robots(i,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        text(robots(i,1)+0.2, robots(i,2), num2str(i));
    end
    plot(formation(:,1), formation(:,2), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
    axis([-5 20 -10 10])
    pause(0.1);
end

%% Cost Function
function J = cost_function(u, x0, y0, theta0, target, T, dt)
    J = 0;
    x = x0; y = y0; theta = theta0;
    lambda = 0.1; % Control effort penalty

    for k = 1:T
        v = u(2*k-1);
        w = u(2*k);
        
        % Update state
        x = x + v * cos(theta) * dt;
        y = y + v * sin(theta) * dt;
        theta = theta + w * dt;
        
        % Cost function: distance to target + control effort
        J = J + norm([x, y] - target)^2 + lambda * (v^2 + w^2);
    end
end

%% Collision Avoidance Constraints
function [c, ceq] = collision_constraints(u, robots, i, T, dt, d_min)
    c = [];
    ceq = [];
    
    x = robots(i, 1);
    y = robots(i, 2);
    theta = robots(i, 3);
    
    % Predict trajectory
    for k = 1:T
        v = u(2*k-1);
        w = u(2*k);
        
        x = x + v * cos(theta) * dt;
        y = y + v * sin(theta) * dt;
        theta = theta + w * dt;
        
        % Check against all other robots
        for j = 1:size(robots,1)
            if j ~= i
                x_other = robots(j,1);
                y_other = robots(j,2);
                
                % Collision constraint: Distance should be greater than d_min
                c = [c; d_min - norm([x - x_other, y - y_other])];
            end
        end
    end
end
