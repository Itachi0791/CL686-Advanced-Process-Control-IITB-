clc; clear; close all;

%% Simulation Parameters
dt = 0.1;            % Time step
T = 10;              % Prediction horizon
N = 5;               % Number of robots
sim_time = 500;      % Total simulation steps
v_max = 5.0;         % Max linear velocity
w_max = pi/2;        % Max angular velocity
d_min = 0.5;         % Minimum safe distance

% Initial positions [x, y, theta]
robots = [0, 0, 0;
          1, 1, pi/4;
         -1, 1, -pi/4;
          1, -1, pi/4;
         -1, -1, -pi/4];

% Formation offsets (relative positions)
formation_offsets = [0, 0;    
                     1, 0.5;  
                    -1, 0.5;
                     1, -0.5;
                    -1, -0.5];

% Store trajectories
trajectories = cell(N,1);
for i = 1:N
    trajectories{i} = robots(i,1:2);
end

%% Coupled MPC Optimization
for t = 1:sim_time
    % Define leader's star trajectory
    theta_t = 0.01 * t; % Parametric time
    leader_x = 5 + 2 * cos(2 * theta_t) * cos(theta_t);
    leader_y = 5 + 2 * cos(2 * theta_t) * sin(theta_t);
    
    % Compute dynamic formation positions
    formation = [leader_x, leader_y;
                 leader_x + formation_offsets(2,1), leader_y + formation_offsets(2,2);
                 leader_x + formation_offsets(3,1), leader_y + formation_offsets(3,2);
                 leader_x + formation_offsets(4,1), leader_y + formation_offsets(4,2);
                 leader_x + formation_offsets(5,1), leader_y + formation_offsets(5,2)];
    
    % Set optimization variables: [v1_1, w1_1, ..., vN_1, wN_1, ..., v1_T, w1_T, ..., vN_T, wN_T]
    u0 = zeros(2*T*N, 1);
    
    % Cost function
    cost_fun = @(u) coupled_cost_function(u, robots, formation, T, dt, N);
    
    % Constraints
    A = []; b = [];
    Aeq = []; beq = [];
    lb = repmat([-v_max; -w_max], T*N, 1);
    ub = repmat([v_max; w_max], T*N, 1);
    
    % Solve coupled MPC optimization
    options = optimoptions('fmincon', 'Display', 'none');
    u_opt = fmincon(cost_fun, u0, A, b, Aeq, beq, lb, ub, ...
                    @(u) coupled_constraints(u, robots, T, dt, d_min, N), options);
    
    % Apply first control input for each robot
    for i = 1:N
        v = u_opt(2*T*(i-1) + 1);
        w = u_opt(2*T*(i-1) + 2);
        robots(i, 1) = robots(i, 1) + v * cos(robots(i, 3)) * dt;
        robots(i, 2) = robots(i, 2) + v * sin(robots(i, 3)) * dt;
        robots(i, 3) = robots(i, 3) + w * dt;
        
        % Store trajectory
        trajectories{i} = [trajectories{i}; robots(i,1:2)];
    end

    % Plotting
    clf; hold on; axis([0 10 0 10]); grid on;
    for i = 1:N
        plot(trajectories{i}(:,1), trajectories{i}(:,2), 'b');
        plot(robots(i,1), robots(i,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        text(robots(i,1)+0.2, robots(i,2), num2str(i));
    end
    plot(formation(:,1), formation(:,2), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
    pause(0.01);
end

%% Coupled Cost Function
function J = coupled_cost_function(u, robots, formation, T, dt, N)
    J = 0;
    lambda = 0.1; % Control effort penalty
    
    for i = 1:N
        x = robots(i, 1);
        y = robots(i, 2);
        theta = robots(i, 3);
        
        for k = 1:T
            v = u(2*T*(i-1) + 2*k-1);
            w = u(2*T*(i-1) + 2*k);
            
            % Predict state
            x = x + v * cos(theta) * dt;
            y = y + v * sin(theta) * dt;
            theta = theta + w * dt;
            
            % Formation constraint: stay close to target position
            J = J + norm([x, y] - formation(i,:))^2 + lambda * (v^2 + w^2);
        end
    end
end

%% Coupled Collision Constraints
function [c, ceq] = coupled_constraints(u, robots, T, dt, d_min, N)
    c = [];
    ceq = [];
    
    % Predict positions for all robots
    predicted_positions = zeros(N, 2);
    for i = 1:N
        x = robots(i, 1);
        y = robots(i, 2);
        theta = robots(i, 3);
        
        for k = 1:T
            v = u(2*T*(i-1) + 2*k-1);
            w = u(2*T*(i-1) + 2*k);
            
            x = x + v * cos(theta) * dt;
            y = y + v * sin(theta) * dt;
            theta = theta + w * dt;
        end
        predicted_positions(i,:) = [x, y];
    end
    
    % Collision avoidance constraints
    for i = 1:N
        for j = i+1:N
            c = [c; d_min - norm(predicted_positions(i,:) - predicted_positions(j,:))];
        end
    end
end
