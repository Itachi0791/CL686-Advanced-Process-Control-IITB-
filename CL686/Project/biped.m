clc; clear; close all;

% Parameters
n = 10;          % Number of states (5 angles + 5 angular velocities)
m = 5;           % Number of control inputs (5 joint torques)
Ns = 200;        % Number of time steps
T = 0.05;        % Sampling time (s)
p = 20;          % Prediction horizon
q = 20;          % Control horizon
tspan = 0:T:Ns*T;

% Define params structure
params.g = 9.81; % Gravity (m/s^2)
params.L_torso = 1; % Torso length (m)
params.L_thigh = 0.5; % Thigh length (m)
params.L_shin = 0.5; % Shin length (m)
params.mass = 1; % Mass (kg)

% Initial state (angles and angular velocities)
x0 = [0; 0.1; -0.2; 0.1; -0.2; zeros(5, 1)]; % Small initial angles, zero velocities

% Weights and constraints
wx = eye(n);    % State weights matrix
wu = 0.001 * eye(m); % Control weights matrix
U_L = -50 * ones(m, 1); U_H = 50 * ones(m, 1); % Control input constraints
delU_L = -10 * ones(m, 1); delU_H = 10 * ones(m, 1); % Control change constraints
X_L = [-pi; -pi/2; -pi/2; -pi/2; -pi/2; -10; -10; -10; -10; -10]; % State constraints (lower)
X_H = [pi; pi/2; pi/2; pi/2; pi/2; 10; 10; 10; 10; 10]; % State constraints (upper)
Xs = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % Setpoint (upright position)

%% MPC Constant matrices

Wx = kron(eye(p), wx); % State weights matrix
Wu = kron(eye(p), wu); % Control weights matrix
Lambda = zeros(m * p, m * p); Lambda_0 = [eye(m); zeros((p - 1) * m, m)]; % Matrices in change of control input inequalities
Lambda(1:m, 1:m) = eye(m);
for k = 2:p
    Lambda((k - 1) * m + 1:k * m, (k - 2) * m + 1:(k - 1) * m) = -eye(m);
    Lambda((k - 1) * m + 1:k * m, (k - 1) * m + 1:k * m) = eye(m);
end
Utilde_L = repmat(U_L, p, 1); Utilde_H = repmat(U_H, p, 1); % Control constraints stacked
Xtilde_L = repmat(X_L, p, 1); Xtilde_H = repmat(X_H, p, 1); % State constraints stacked
DeltaU_L = repmat(delU_L, p, 1); DeltaU_H = repmat(delU_H, p, 1); % Control change constraints stacked
Zs = repmat(Xs, p, 1);
Aeq = zeros(m * (p - q), m * p); 
for i = 1:(p - q)
    Aeq((i - 1) * m + 1:i * m, (q - 1) * m + 1:q * m) = eye(m);
    Aeq((i - 1) * m + 1:i * m, (q + i - 1) * m + 1:(q + i) * m) = -eye(m);
end    
beq = zeros(m * (p - q), 1);
A_inequality = @(Su) [eye(m * p); Lambda; Su; -eye(m * p); -Lambda; -Su];
b_inequality = @(x_k, u_kminus1, Del_k, Sx, Sdel) [Utilde_H; DeltaU_H + Lambda_0 * u_kminus1; Xtilde_H - Sx * x_k - Sdel * Del_k;...
                                         -Utilde_L; -DeltaU_L - Lambda_0 * u_kminus1; -Xtilde_L + Sx * x_k + Sdel * Del_k];

%% MPC simulation

X = zeros(n, Ns + 1); X(:, 1) = x0; U = zeros(m, length(tspan) - 1);
Ufk_prev = zeros(m * p, 1);
options = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off');
fprintf("Starting MPC Computation : \n")
for k = 2:Ns
    % To display progress
    if mod(k, round(Ns / 5)) == 0
        fprintf('%d / %d iterations done\n', k, Ns);
    end
    B_mat = B_matrix(X(:, k), params);
    A_mat = A_matrix(X(:, k), params);
    Phi_mat = Phi_matrix(A_mat, T);
    Gamma_mat = Gamma_matrix(A_mat, B_mat, T, n, m);
    del_k = T * biped_dyn(X(:, k), params, U(:, k - 1)) + X(:, k) - Phi_mat * X(:, k) - Gamma_mat * U(:, k - 1);
    Del_k = repmat(del_k, p, 1);
    [Sx, Su, Sdel, H, F] = MPC_matrices(Phi_mat, Gamma_mat, p, Wx, Wu, X(:, k), Del_k, Zs);
    A_ineq = A_inequality(Su);
    b_ineq = b_inequality(X(:, k), U(:, k - 1), Del_k, Sx, Sdel);
    Ufk = quadprog(H, F, A_ineq, b_ineq, Aeq, beq, [], [], Ufk_prev, options);
    U(:, k) = Lambda_0' * Ufk;
    [~, Y] = ode45(@(t, x) biped_dyn(x, params, U(:, k)), [0, T], X(:, k)); 
    X(:, k + 1) = (Y(end, :))';
    Ufk_prev = Ufk;
end
fprintf("MPC Computation completed.\n")

%% Animation

figure;
hold on;
axis equal;
grid on;
xlim([-2, 2]);
ylim([-2, 2]);
title('Biped Robot Animation');

% Link lengths
L_torso = params.L_torso;
L_thigh = params.L_thigh;
L_shin = params.L_shin;

for k = 1:Ns
    % Extract joint angles
    q1 = X(1, k); % Torso
    q2 = X(2, k); % Left thigh
    q3 = X(3, k); % Left shin
    q4 = X(4, k); % Right thigh
    q5 = X(5, k); % Right shin
    
    % Compute positions
    torso_x = [0, L_torso * sin(q1)];
    torso_y = [0, L_torso * cos(q1)];
    
    left_thigh_x = [torso_x(2), torso_x(2) + L_thigh * sin(q1 + q2)];
    left_thigh_y = [torso_y(2), torso_y(2) + L_thigh * cos(q1 + q2)];
    
    left_shin_x = [left_thigh_x(2), left_thigh_x(2) + L_shin * sin(q1 + q2 + q3)];
    left_shin_y = [left_thigh_y(2), left_thigh_y(2) + L_shin * cos(q1 + q2 + q3)];
    
    right_thigh_x = [torso_x(2), torso_x(2) + L_thigh * sin(q1 + q4)];
    right_thigh_y = [torso_y(2), torso_y(2) + L_thigh * cos(q1 + q4)];
    
    right_shin_x = [right_thigh_x(2), right_thigh_x(2) + L_shin * sin(q1 + q4 + q5)];
    right_shin_y = [right_thigh_y(2), right_thigh_y(2) + L_shin * cos(q1 + q4 + q5)];
    
    % Plot
    cla;
    plot(torso_x, torso_y, 'k', 'LineWidth', 3); % Torso
    plot(left_thigh_x, left_thigh_y, 'b', 'LineWidth', 2); % Left thigh
    plot(left_shin_x, left_shin_y, 'b', 'LineWidth', 2); % Left shin
    plot(right_thigh_x, right_thigh_y, 'r', 'LineWidth', 2); % Right thigh
    plot(right_shin_x, right_shin_y, 'r', 'LineWidth', 2); % Right shin
    drawnow;
    pause(0.05);
end

%% Functions

function A = A_matrix(x, params)
    % Linearized state matrix (approximate)
    A = zeros(10, 10);
    A(1:5, 6:10) = eye(5);
    A(6:10, 1:5) = diag([-params.g, -params.g, -params.g, -params.g, -params.g]); % Gravity effect
end

function B = B_matrix(x, params)
    % Linearized input matrix (approximate)
    B = [zeros(5, 5); eye(5)];
end

function Phi = Phi_matrix(A_matrix, T)
    Phi = expm(A_matrix * T);
end

function Gamma = Gamma_matrix(A_matrix, B, T, n, m)
    Gamma = (expm(A_matrix * T) - eye(n)) * pinv(A_matrix) * B;
end

function xdot = biped_dyn(x, params, u)
    % Nonlinear dynamics of the biped robot
    q = x(1:5);
    qdot = x(6:10);
    qddot = -params.g * sin(q) + u; % Simplified dynamics (gravity + control input)
    xdot = [qdot; qddot];
end