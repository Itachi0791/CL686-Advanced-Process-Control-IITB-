clc; clear; close all;

% CartPole parameters
m_p = 1;    % Mass of the pole (kg)
m_c = 4.0;  % Mass of the cart (kg)
l = 1;      % Length of the pole (m)
g = 9.81;   % Acceleration due to gravity (m/s^2)
b = 0.1;    % Damping coefficient (optional)

params.m_p = m_p;
params.m_c = m_c;
params.l = l;
params.g = g;
params.b = b;

% MPC parameters
n = 4;     % Number of states [x; theta; x_dot; theta_dot]
m = 1;     % Number of inputs [u]
Ns = 200;  % No. of time steps
T = 0.05;  % Sampling Time
p = 100;   % Prediction Horizon   
q = 20;    % Control Horizon
tspan = 0:T:Ns*T;

% Initial state and reference
x0 = [0; pi; 0; 0]; % [x; theta; x_dot; theta_dot] - pendulum starts upside down
wx = diag([10, 50, 1, 1]); % State weights matrix - higher weight on position and angle
wu = 0.001*eye(m);   % Control weights matrix

% Constraints
U_L = -5*ones(m,1); U_H = 5*ones(m,1);          % Control input constraints
delU_L = -1*ones(m,1); delU_H = 1*ones(m,1);      % Control change constraints
X_L = [-2; -pi; -2; -5]; X_H = [2; pi; 2; 5];     % State constraints
Xs = [0; 0; 0; 0];                                % Setpoint - upright position at origin

%% MPC Constant matrices
Wx = kron(eye(p), wx);    % State weights for prediction horizon
Wu = kron(eye(p), wu);    % Control weights for prediction horizon

% Matrices for change of control input inequalities
Lambda = zeros(m*p, m*p); 
Lambda_0 = [eye(m); zeros((p-1)*m, m)];
Lambda(1:m, 1:m) = eye(m);
for k = 2:p
    Lambda((k-1)*m + 1:k*m, (k-2)*m + 1:(k-1)*m) = -eye(m);
    Lambda((k-1)*m + 1:k*m, (k-1)*m + 1:k*m) = eye(m);
end

% Stacked constraints
Utilde_L = repmat(U_L, p, 1); Utilde_H = repmat(U_H, p, 1);          % Control constraints
Xtilde_L = repmat(X_L, p, 1); Xtilde_H = repmat(X_H, p, 1);          % State constraints
DeltaU_L = repmat(delU_L, p, 1); DeltaU_H = repmat(delU_H, p, 1);    % Control change constraints
Zs = repmat(Xs, p, 1);                                               % Stacked setpoint

% Equality constraints for control horizon
Aeq = zeros(m*(p-q), m*p); 
for i = 1:(p-q)
    Aeq((i-1)*m+1:i*m, (q-1)*m+1:q*m) = eye(m);
    Aeq((i-1)*m+1:i*m, (q+i-1)*m+1:(q+i)*m) = -eye(m);
end    
beq = zeros(m*(p-q), 1);

% Inequality constraints functions
A_inequality = @(Su)[eye(m*p); Lambda; Su; -eye(m*p); -Lambda; -Su];
b_inequality = @(x_k, u_kminus1, Del_k, Sx, Sdel) [Utilde_H; DeltaU_H+Lambda_0*u_kminus1; Xtilde_H-Sx*x_k-Sdel*Del_k;...
                                      -Utilde_L; -DeltaU_L-Lambda_0*u_kminus1; -Xtilde_L+Sx*x_k+Sdel*Del_k];

%% MPC simulation
X = zeros(n, Ns+1); X(:,1) = x0; 
U = zeros(m, length(tspan)-1);

% First step simulation with U(0) = 0
[~, Y] = ode45(@(t,x) cartpole_dyn(t, x, params, U(:,1)), [0, T], X(:,1)); 
X(:,2) = (Y(end,:))';

Ufk_prev = zeros(m*p, 1);
options = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off');

fprintf("Starting MPC Computation : \n")
for k = 2:Ns
    % To display progress
    if mod(k, round(Ns/5)) == 0
        fprintf('%d / %d iterations done\n', k, Ns);
    end
    
    % Linearize the system at current state and input
    [A_mat, B_mat] = linearize_cartpole(X(:,k), U(:,k-1), params);
    
    % Discretize the system
    Phi_mat = Phi_matrix(A_mat, T);
    Gamma_mat = Gamma_matrix(A_mat, B_mat, T, n, m);
    
    % Calculate the disturbance term
    x_dot_nonlinear = cartpole_dyn(0, X(:,k), params, U(:,k-1));
    x_dot_linear = A_mat * X(:,k) + B_mat * U(:,k-1);
    del_k = T * (x_dot_nonlinear - x_dot_linear);
    Del_k = repmat(del_k, p, 1);
    
    % Calculate MPC matrices
    [Sx, Su, Sdel, H, F] = MPC_matrices(Phi_mat, Gamma_mat, p, Wx, Wu, X(:,k), Del_k, Zs);
    
    % Solve quadratic program
    A_ineq = A_inequality(Su);
    b_ineq = b_inequality(X(:,k), U(:,k-1), Del_k, Sx, Sdel);
    Ufk = quadprog(H, F, A_ineq, b_ineq, Aeq, beq, [], [], Ufk_prev, options);
    
    % Apply first control input
    U(:,k) = Lambda_0' * Ufk;
    
    % Simulate system for one step
    [~, Y] = ode45(@(t,x) cartpole_dyn(t, x, params, U(:,k)), [0, T], X(:,k)); 
    X(:,k+1) = (Y(end,:))';
    
    % Update previous control sequence
    Ufk_prev = Ufk;
end
fprintf("MPC Computation completed.\n")

%% Plots
figure;
sgtitle("Adaptive MPC - Controlled States", "FontSize", 25, "FontWeight", "bold");

% Plot cart position
subplot(2,2,1)
plot(tspan, X(1,:), 'LineWidth', 2);
hold on;
plot(tspan, Xs(1)*ones(size(tspan)), "LineStyle", "-.", 'LineWidth', 2);
plot(tspan, X_H(1)*ones(size(tspan)), "k", "LineStyle", "--", 'LineWidth', 2);
plot(tspan, X_L(1)*ones(size(tspan)), "k", "LineStyle", "--", 'LineWidth', 2);
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold')
legend('Actual', 'Reference', 'Constraints', 'FontSize', 12)
grid on
set(gca, 'fontsize', 14)

% Plot pendulum angle
subplot(2,2,2)
plot(tspan, X(2,:), 'LineWidth', 2);
hold on;
plot(tspan, Xs(2)*ones(size(tspan)), "LineStyle", "-.", 'LineWidth', 2);
plot(tspan, X_H(2)*ones(size(tspan)), "k", "LineStyle", "--", 'LineWidth', 2);
plot(tspan, X_L(2)*ones(size(tspan)), "k", "LineStyle", "--", 'LineWidth', 2);
ylabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold')
legend('Actual', 'Reference', 'Constraints', 'FontSize', 12)
grid on
set(gca, 'fontsize', 14)

% Plot cart velocity
subplot(2,2,3)
plot(tspan, X(3,:), 'LineWidth', 2);
hold on;
plot(tspan, Xs(3)*ones(size(tspan)), "LineStyle", "-.", 'LineWidth', 2);
plot(tspan, X_H(3)*ones(size(tspan)), "k", "LineStyle", "--", 'LineWidth', 2);
plot(tspan, X_L(3)*ones(size(tspan)), "k", "LineStyle", "--", 'LineWidth', 2);
ylabel('$\dot{x}$', 'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold')
legend('Actual', 'Reference', 'Constraints', 'FontSize', 12)
grid on
set(gca, 'fontsize', 14)

% Plot pendulum angular velocity
subplot(2,2,4)
plot(tspan, X(4,:), 'LineWidth', 2);
hold on;
plot(tspan, Xs(4)*ones(size(tspan)), "LineStyle", "-.", 'LineWidth', 2);
plot(tspan, X_H(4)*ones(size(tspan)), "k", "LineStyle", "--", 'LineWidth', 2);
plot(tspan, X_L(4)*ones(size(tspan)), "k", "LineStyle", "--", 'LineWidth', 2);
xlabel("Time [s]", 'FontSize', 17)
ylabel('$\dot{\theta}$', 'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold')
legend('Actual', 'Reference', 'Constraints', 'FontSize', 12)
grid on
set(gca, 'fontsize', 14)

% Plot control input
figure;
stairs(tspan(1:end-1), U, 'LineWidth', 2)
hold on
plot(tspan, U_H(1)*ones(size(tspan)), "k", "LineStyle", "--", 'LineWidth', 2);
plot(tspan, U_L(1)*ones(size(tspan)), "k", "LineStyle", "--", 'LineWidth', 2);
ylabel('Control Input', 'FontSize', 20, 'FontWeight', 'bold')
xlabel("Time [s]", 'FontSize', 17)
legend('Actual', 'Constraints', 'FontSize', 12)
title("Adaptive MPC - Manipulated Inputs", "FontSize", 25, "FontWeight", "bold");
grid on
set(gca, 'fontsize', 18)

%% Animate cartpole system
figure('Position', [100, 100, 800, 600]);
hold on;

% Initialize cartpole graphics
cart_width = 0.4;
cart_height = 0.2;
cart = rectangle('Position', [-cart_width/2, -cart_height/2, cart_width, cart_height], 'FaceColor', 'b');
pole = line([0, 0], [0, 0], 'Color', 'k', 'LineWidth', 2);
bob = plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
time_text = text(-2, 2, 'Time: 0.00 s', 'FontSize', 14, 'Color', 'r');

% Set axis limits and aspect ratio
axis equal;
grid on;
xlim([-2.5, 2.5]);
ylim([-1, 2]);
xlabel('Position (m)', 'FontSize', 14);
ylabel('Height (m)', 'FontSize', 14);
title('CartPole with Adaptive MPC Control', 'FontSize', 16);
set(gca, 'FontSize', 12);

% Animation loop
for k = 1:10:length(tspan)
    % Extract state variables
    x_pos = X(1, k);
    theta = X(2, k);
    
    % Calculate pendulum end position
    pendulum_x = x_pos + l * sin(theta);
    pendulum_y = l * cos(theta);
    
    % Update graphics
    set(cart, 'Position', [x_pos-cart_width/2, -cart_height/2, cart_width, cart_height]);
    set(pole, 'XData', [x_pos, pendulum_x], 'YData', [0, pendulum_y]);
    set(bob, 'XData', pendulum_x, 'YData', pendulum_y);
    set(time_text, 'String', sprintf('Time: %.2f s', (k-1)*T));
    
    % Pause to control animation speed
    drawnow;
    pause(0.01);
end

hold off;

%% Helper Functions

function [A, B] = linearize_cartpole(x, u, params)
    % Linearize the cartpole dynamics at the current state and input
    % Returns state-space matrices A and B for linearized system
    
    % Extract parameters
    m_p = params.m_p;
    m_c = params.m_c;
    l = params.l;
    g = params.g;
    b = params.b;
    
    % Extract states
    theta = x(2);
    theta_dot = x(4);
    
    % Precompute trig functions
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    
    % Define the mass ratio
    m_total = m_c + m_p;
    
    % Define common denominators
    den = l * (4/3 - m_p * (cos_theta^2) / m_total);
    
    % Initialize matrices
    A = zeros(4, 4);
    B = zeros(4, 1);
    
    % Fill in A matrix
    % First row: dx/dt = x_dot
    A(1, 3) = 1;
    
    % Second row: dtheta/dt = theta_dot
    A(2, 4) = 1;
    
    % Third row: dx_dot/dt partial derivatives (cart acceleration)
    A(3, 2) = m_p * l * (theta_dot^2 * cos_theta + g * sin_theta * cos_theta / l) / m_total;
    if abs(theta) < 1e-3  % Near the upright position, use linearized form
        A(3, 2) = m_p * g / m_total;
    end
    A(3, 3) = -b / m_total;
    A(3, 4) = -m_p * l * sin_theta * theta_dot / m_total;
    
    % Fourth row: dtheta_dot/dt partial derivatives (angular acceleration)
    A(4, 2) = g * cos_theta / den + m_p * cos_theta * sin_theta / (den * m_total);
    A(4, 3) = -b * cos_theta / (m_total * den);
    A(4, 4) = -m_p * l * cos_theta * sin_theta * theta_dot / (den * m_total);
    
    % Fill in B matrix (control input effects)
    B(3, 1) = 1 / m_total;
    B(4, 1) = -cos_theta / (m_total * den);

end

function x_dot = cartpole_dyn(t, x, params, u)
    % Full nonlinear dynamics of cartpole system
    % State x = [x; theta; x_dot; theta_dot]
    % u = force applied to cart
    
    % Extract parameters
    m_p = params.m_p;
    m_c = params.m_c;
    l = params.l;
    g = params.g;
    b = params.b;
    
    % Extract states
    theta = x(2);
    x_dot = x(3);
    theta_dot = x(4);
    
    % Compute trig values
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    
    % Compute common denominators
    den1 = m_c + m_p - m_p * cos_theta^2;
    den2 = l * (4/3 - m_p * cos_theta^2 / (m_c + m_p));
    
    % Compute accelerations
    x_ddot = (u + m_p * sin_theta * (l * theta_dot^2 + g * cos_theta) - b * x_dot) / den1;
    theta_ddot = (-u * cos_theta - m_p * l * theta_dot^2 * sin_theta * cos_theta - (m_c + m_p) * g * sin_theta - b * cos_theta * x_dot) / ((m_c + m_p) * den2);
    
    % Return state derivative
    x_dot = [x_dot; theta_dot; x_ddot; theta_ddot];
end

function Phi = Phi_matrix(A_matrix, T)
    % Compute the state transition matrix
    Phi = expm(A_matrix * T);
end

function Gamma = Gamma_matrix(A_matrix, B, T, n, m)
    % Compute the input effect matrix
    if det(A_matrix) > 1e-10
        Gamma = (expm(A_matrix * T) - eye(n)) * inv(A_matrix) * B;
    else
        % Use numerical integration for singular A
        ode_func = @(tau, x) reshape(expm(A_matrix * tau) * B, [], 1);
        [~, Y] = ode45(ode_func, [0, T], zeros(n * m, 1));
        Gamma = reshape(Y(end, :)', n, m);
    end
end
