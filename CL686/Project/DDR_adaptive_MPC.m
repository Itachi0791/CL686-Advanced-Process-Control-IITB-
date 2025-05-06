clc; clear; close all;

%% Parameters
T = 0.1; % Sampling time
Ns = 500; % Number of time steps
p = 80; % Prediction horizon
q = 50; % Control horizon
n = 3; m =2; % State and input dimensions
wx = diag([10, 10, 0]); % State weights
wu = diag([0.1, 0.1]); % Control weights
U_L = [-5; -5]; % Lower control input constraints
U_H = [5; 5]; % Upper control input constraints
delU_L = [-1; -1]; % Lower control change constraints
delU_H = [1; 1]; % Upper control change constraints
X_L = [-10; -10; -pi]; % Lower state constraints
X_H = [10; 10; pi]; % Upper state constraints

% Circular trajectory
t = linspace(0, 2*pi, Ns+1);
Xs = [-5*ones(size(t));6*ones(size(t));zeros(size(t))];%[4*cos(t); 8*sin(t); zeros(1, Ns+1)]; % Setpoint (circular trajectory)

% Initial state
x0 = [0;0;atan2(-5,6);];% Xs(:,1);

%% MPC Constant matrices

Wx = kron(eye(p),wx); % using kronecker product we can stack wx in block diagonal matrix p times
Wu = kron(eye(p),wu); % same for wu
Lambda = zeros(m*p,m*p); Lambda_0 = [eye(m);zeros((p-1)*m,m)]; % Matrices in change of control input inequalities
Lambda(1:m,1:m) = eye(m);
for k = 2:p
    Lambda((k-1)*m + 1:k*m,(k-2)*m + 1:(k-1)*m) = -eye(m);
    Lambda((k-1)*m + 1:k*m,(k-1)*m + 1:k*m) = eye(m);
end
Utilde_L = repmat(U_L,p,1); Utilde_H = repmat(U_H,p,1); % Control constraints stacked
Xtilde_L = repmat(X_L,p,1); Xtilde_H = repmat(X_H,p,1); % State constraints stacked
DeltaU_L = repmat(delU_L,p,1); DeltaU_H = repmat(delU_H,p,1); % Control change constraints stacked
Zs = repmat(Xs,p,1);
Aeq = zeros(m*(p-q),m*p); 
for i = 1:(p-q)
    Aeq((i-1)*m+1:i*m,(q-1)*m+1:q*m) = eye(m);
    Aeq((i-1)*m+1:i*m,(q+i-1)*m+1:(q+i)*m) = -eye(m);
end    
beq = zeros(m*(p-q),1);
A_inequality = @(Su)[eye(m*p);Lambda;Su;-eye(m*p);-Lambda;-Su];
b_inequality = @(x_k,u_kminus1,Del_k,Sx,Sdel) [Utilde_H;DeltaU_H+Lambda_0*u_kminus1;Xtilde_H-Sx*x_k-Sdel*Del_k;...
                                         -Utilde_L;-DeltaU_L-Lambda_0*u_kminus1;-Xtilde_L+Sx*x_k+Sdel*Del_k];

%% MPC simulation
tspan = 0:T:Ns*T;
X = zeros(n,Ns+1); X(:,1)=x0; U = zeros(m,length(tspan)-1);
% First step simulation with U(0) = 0
[~,Y] = ode45(@(t,x) diff_drive_dyn(x,U(:,1)),[0,T],X(:,1)); 
X(:,2) = (Y(end,:))';
Ufk_prev = zeros(m*p,1);
options = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off');
fprintf("Starting MPC Computation : \n")
for k = 2:Ns
    if mod(k, round(Ns/5)) == 0
        fprintf('%d / %d iterations done\n', k, Ns);
    end
    A_mat = A_matrix(X(:, k));
    B_mat = B_matrix(X(:, k));
    Phi_mat = Phi_matrix(A_mat, T);
    Gamma_mat = Gamma_matrix(A_mat, B_mat, T, 3, 2);
    del_k = T * diff_drive_dyn(X(:, k), U(:, k-1)) + X(:, k) - Phi_mat * X(:, k) - Gamma_mat * U(:, k-1);
    Del_k = repmat(del_k, p, 1);
    Zs = repmat(Xs(:,k), p, 1);
    [Sx, Su, Sdel, H, F] = MPC_matrices(Phi_mat, Gamma_mat, p, Wx, Wu, X(:, k), Del_k, Zs);
    A_ineq = A_inequality(Su);
    b_ineq = b_inequality(X(:, k), U(:, k-1), Del_k, Sx, Sdel);
    Ufk = quadprog(H, F, A_ineq, b_ineq, Aeq, beq, [], [], Ufk_prev, options);
    U(:, k) = Lambda_0' * Ufk;
    [~,Y] = ode45(@(t,x) diff_drive_dyn(x,U(:,k)),[0,T],X(:,k)); 
    X(:,k+1) = (Y(end,:))';
    Ufk_prev = Ufk;
end
fprintf("MPC Computation completed.\n")

%% Plots
figure;
sgtitle("Adaptive MPC - Controlled States", "FontSize", 25, "FontWeight", "bold");
subplot(3, 1, 1)
plot(0:T:Ns*T, X(1, :), 'LineWidth', 2);
hold on;
plot(0:T:Ns*T, Xs(1, :), 'LineStyle', '-.', 'LineWidth', 2);
plot(0:T:Ns*T, X_H(1)*ones(1, Ns+1), 'k', 'LineStyle', '--', 'LineWidth', 2);
plot(0:T:Ns*T, X_L(1)*ones(1, Ns+1), 'k', 'LineStyle', '--', 'LineWidth', 2);
ylabel('x', 'FontSize', 25, 'FontWeight', 'bold', 'Rotation', 0)
legend('Actual', 'Reference', 'Constraints', 'FontSize', 12)
grid on
set(gca, 'fontsize', 18)
subplot(3, 1, 2)
plot(0:T:Ns*T, X(2, :), 'LineWidth', 2);
hold on;
plot(0:T:Ns*T, Xs(2, :), 'LineStyle', '-.', 'LineWidth', 2);
plot(0:T:Ns*T, X_H(2)*ones(1, Ns+1), 'k', 'LineStyle', '--', 'LineWidth', 2);
plot(0:T:Ns*T, X_L(2)*ones(1, Ns+1), 'k', 'LineStyle', '--', 'LineWidth', 2);
ylabel('y', 'FontSize', 25, 'FontWeight', 'bold', 'Rotation', 0)
legend('Actual', 'Reference', 'Constraints', 'FontSize', 12)
grid on
set(gca, 'fontsize', 18)
subplot(3, 1, 3)
plot(0:T:Ns*T, X(3, :), 'LineWidth', 2);
hold on;
plot(0:T:Ns*T, Xs(3, :), 'LineStyle', '-.', 'LineWidth', 2);
plot(0:T:Ns*T, X_H(3)*ones(1, Ns+1), 'k', 'LineStyle', '--', 'LineWidth', 2);
plot(0:T:Ns*T, X_L(3)*ones(1, Ns+1), 'k', 'LineStyle', '--', 'LineWidth', 2);
xlabel("Time [in s]", 'FontSize', 17)
ylabel('\theta', 'FontSize', 25, 'FontWeight', 'bold', 'Rotation', 0)
legend('Actual', 'Reference', 'Constraints', 'FontSize', 12)
grid on
set(gca, 'fontsize', 18)

figure;
stairs(0:T:(Ns-1)*T, U', 'LineWidth', 2)
hold on
plot(0:T:Ns*T, U_H(1)*ones(1, Ns+1), 'k', 'LineStyle', '--', 'LineWidth', 2);
plot(0:T:Ns*T, U_L(1)*ones(1, Ns+1), 'k', 'LineStyle', '--', 'LineWidth', 2);
ylabel('U', 'FontSize', 25, 'FontWeight', 'bold', 'Rotation', 0)
xlabel("Time [in s]", 'FontSize', 17)
legend('Actual', 'Constraints', 'FontSize', 12)
title("Adaptive MPC - Manipulated Inputs", "FontSize", 25, "FontWeight", "bold");
grid on
set(gca, 'fontsize', 18)

% Animation
figure('Position', [100, 100, 700, 700]);
hold on;
axis equal;
grid on;
axis([-10, 10, -10, 10]);
xlabel('X axis (m)');
ylabel('Y axis (m)');
title('Adaptive MPC for Differential Drive Robot');
plot(Xs(1,:), Xs(2,:), 'g--', 'LineWidth', 1.5); % Reference trajectory

% Robot parameters
wheel_radius = 0.1;
wheel_width = 0.1;
car_width = 0.6;
car_length = 1.0;

% Car body shape (before rotation)
car_vertices = [-car_length/2, -car_width/2;
                car_length/2, -car_width/2;
                car_length/2, car_width/2;
               -car_length/2, car_width/2]';

% Wheel shape (before rotation)
wheel_vertices = [-wheel_radius/2, -wheel_width/2;
                  wheel_radius/2, -wheel_width/2;
                  wheel_radius/2, wheel_width/2;
                 -wheel_radius/2, wheel_width/2]';

% Initialize car body patch
car_body_patch = patch('XData', [], 'YData', [], 'FaceColor', 'b');

% Initialize wheels
front_left_wheel_patch = patch('XData', [], 'YData', [], 'FaceColor', 'k');
front_right_wheel_patch = patch('XData', [], 'YData', [], 'FaceColor', 'k');
rear_left_wheel_patch = patch('XData', [], 'YData', [], 'FaceColor', 'k');
rear_right_wheel_patch = patch('XData', [], 'YData', [], 'FaceColor', 'k');

% Initialize time text
time_text = text(-1.8, 1.8, 'Time: 0.00 s', 'FontSize', 18, 'Color', 'r');

% Initialize robot path
robot_path = plot(X(1,1), X(2,1), 'k:', 'LineWidth', 1.5);

pause(1); % Pause before starting animation

for k = 1:Ns+1
    % Car position and orientation
    x_pos = X(1, k);
    y_pos = X(2, k);
    theta = X(3, k); % Orientation angle

    % Rotation matrix
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

    % Rotate and update car body
    rotated_car = R * car_vertices + [x_pos; y_pos];
    set(car_body_patch, 'XData', rotated_car(1,:), 'YData', rotated_car(2,:));

    % Define wheel center positions
    wheel_positions = [ -car_length/2, car_width/2 - wheel_width;  % Front left
                         car_length/2, car_width/2 - wheel_width;  % Front right
                        -car_length/2, -car_width/2;               % Rear left
                         car_length/2, -car_width/2 ];             % Rear right

    % Rotate and update wheels
    wheel_patches = {front_left_wheel_patch, front_right_wheel_patch, rear_left_wheel_patch, rear_right_wheel_patch};
    for w = 1:4
        rotated_wheel = R * (wheel_vertices + wheel_positions(w,:)') + [x_pos; y_pos];
        set(wheel_patches{w}, 'XData', rotated_wheel(1,:), 'YData', rotated_wheel(2,:));
    end

    % Update time text
    set(time_text, 'String', sprintf('Time: %.2f s', (k-1)*T));

    % Update robot path
    set(robot_path, 'XData', X(1,1:k), 'YData', X(2,1:k));

    % Pause for animation speed
    pause(0.05);
end
hold off;


%% Functions

function A = A_matrix(x)
    A = [0, 0, -sin(x(3)); 0, 0, cos(x(3)); 0, 0, 0];
end

function B = B_matrix(x)
    B = [cos(x(3)), 0; sin(x(3)), 0; 0, 1];
end

function Phi = Phi_matrix(A_matrix, T)
    Phi = expm(A_matrix * T);
end

function Gamma = Gamma_matrix(A_matrix, B, T, n, m)
    if det(A_matrix) > 1e-10
        Gamma = (expm(A_matrix * T) - eye(n)) * pinv(A_matrix) * B;
    else
        ode_func = @(tau, x) reshape(expm(A_matrix * tau) * B, [], 1);
        [~, Y] = ode45(ode_func, [0, T], zeros(n * m, 1));
        Gamma = reshape(Y(end, :)', n, m);
    end
end

function xdot = diff_drive_dyn(x, u)
    xdot = [cos(x(3)), 0; sin(x(3)), 0; 0, 1] * u;
end

