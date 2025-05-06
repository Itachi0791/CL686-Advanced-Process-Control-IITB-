clc; clear; close all;

%% Parameters
dt = 0.1;       % Time step [s]
N = 10;         % Prediction horizon
T_final = 5;    % Total simulation time [s]
g = 9.81;       % Gravity [m/s^2]
h0 = 0.9;       % Initial CoM height [m]

% Initial state: [x, y, vx, vy] (CoM position & velocity)
s = [0; 0; 0.2; 0];  % Start moving in x-direction

% Footstep locations (initial step at (0,0))
footsteps = [0, 0];  

% Optimization weights
Q = diag([1, 1, 0.1, 0.1]);   % State cost
R = diag([0.01, 0.01]);       % Control cost

% Constraints
step_length_max = 0.3;  % Maximum step length [m]
step_width_max = 0.2;   % Maximum step width [m]

% Time vector
t_vec = 0:dt:T_final;

%% Simulation
for t = t_vec
    
    % Compute variable height h_xy based on CoM motion
    h_xy = h0 + 0.05 * sin(2 * pi * t / T_final);  % Example variation
    
    % Linearized Discrete Dynamics (VH-LIPM)
    A_d = [1, 0, dt, 0;
           0, 1, 0, dt;
           dt * g / h_xy, 0, 1, 0;
           0, dt * g / h_xy, 0, 1];
       
    B_d = [0, 0;
           0, 0;
           dt, 0;
           0, dt];  % Control input affects velocity

    % MPC Optimization (Quadprog)
    H = blkdiag(kron(eye(N), Q), kron(eye(N-1), R));  % Quadratic cost matrix
    f = zeros(size(H,1),1);  % Linear cost (zero)

    % Constraint Matrices
    A_cons = kron(eye(N-1), [1, -1]);  % Step constraints
    b_cons = repmat([step_length_max; step_width_max], N-1, 1);

    % Solve Optimization
    options = optimoptions('quadprog', 'Display', 'off');
    u_opt = quadprog(H, f, A_cons, b_cons, [], [], [], [], [], options);
    
    % Extract first step
    u_opt = u_opt(1:2);
    
    % Convert footstep to velocity input
    v_des = (u_opt - s(1:2)) / dt;
    
    % Update state
    s = A_d * s + B_d * v_des;
    
    % Store footstep
    footsteps = [footsteps; u_opt'];

    % Visualization
    clf;
    hold on; grid on; axis equal;
    xlim([-0.5, 2]); ylim([-0.5, 0.5]);  
    plot(footsteps(:,1), footsteps(:,2), 'ro-', 'MarkerSize', 10, 'LineWidth', 2);
    plot(s(1), s(2), 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b'); % CoM
    drawnow;
    
end
