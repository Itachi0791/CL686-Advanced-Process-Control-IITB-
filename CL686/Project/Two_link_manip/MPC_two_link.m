clc; clear; close all;

% Physical parameters
m1 = 5; m2 = 4; l1 = 1; l2 = 0.8; I1 = m1*l1^2/3; I2 = m2*l2^2/3; g = 9.81;
params.m1 = m1; params.m2 = m2; params.l1 = l1;params.l2 = l2; params.I1 =I1; params.I2 = I2; params. g = g;

% Dynamics matrices
H = @(q) [I1 + m1*l1^2/4 + I2 + m2*(l2^2/4 + l1^2 + l1*l2*cos(q(2))), I2 + m2*(l2^2/4 + l1*l2*cos(q(2))/2); 
          I2 + m2*(l2^2/4 + l1*l2*cos(q(2))/2), I2 + m2*l2^2/4];
C = @(q) [-m2*l1*l2*sin(q(2))*q(4), -m2*l1*l2*sin(q(2))*q(4)/2; 
           m2*l1*l2*sin(q(2))*q(3)/2, 0];
G = @(q) [m1*g*l1*cos(q(1))/2 + m2*g*(l1*cos(q(1)) + l2/2*cos(q(1)+q(2))); 
          m2*g*l2/2*cos(q(1)+q(2))];

compute_manipulator_jacobians(m1,m2,l1,l2,g);


% Dynamics function for ode89
arm_dynamics = @(x, u) [x(3:4); H(x) \ (u - C(x)*x(3:4) - G(x))];

%% MPC parameters and constraints

T = 0.05; % Sampling Time (in s)
Ns = 800; % No. of Simulation Steps
tspan = 0:T:Ns*T;
p = 40; % Prediction Horizon   
q = 20; % Control Horizon
n = 4; m = 2; % No. of states and control inputs respectively 
Xs = [3*pi/8*sin(0.4*tspan);pi/8-pi/6*cos(0.3*tspan);3*pi/8*0.4*cos(0.4*tspan);pi/6*0.3*sin(0.3*tspan)]; 
%Xs = [pi*ones(size(tspan));0*pi/2*ones(size(tspan));zeros(size(tspan));zeros(size(tspan))];
Us = find_steady_state_input_trajectory(Xs, m, params);
wx = eye(n); % State Weights matrix
wu = 0.001*eye(m); % Control weights matrix
U_L = -20*ones(m,1); U_H = 20*ones(m,1); % Control input constraints
delU_L = -5*ones(m,1); delU_H = 5*ones(m,1); % Control change constraints
X_L = -100*ones(n,1); X_H = 100*ones(n,1); % State constraints

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
%Zs = repmat(Xs,p,1);
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

X = zeros(n,Ns+1); U = zeros(m,length(tspan)-1);
X0 = Xs(:,1); X(:,1) = X0;% U(:,1) = Us;
options_ode = odeset('abstol', 1e-9, 'reltol', 1e-9);
% First step simulation with U(0) = 0
[~,Y] = ode45(@(t,Z) arm_dynamics(Z,U(:,1)),[0,T],X(:,1),options_ode); 
X(:,2) = (Y(end,:))';
Ufk_prev = zeros(m*p,1);
options = optimoptions('quadprog','Algorithm','active-set','Display','off');
fprintf("Starting MPC Computation : \n")
for k = 2:Ns
    % To display progress
    if mod(k, round(Ns/5)) == 0
        fprintf('%d / %d iterations done\n', k, Ns);
    end
    % Substitute values into A and B
    A_mat = A_matrix(X(:,k),U(:,k-1));
    B_mat = B_matrix(X(:,k),U(:,k-1));
    Phi_mat = expm(A_mat*T);
    Gamma_mat = Gamma_matrix(A_mat, B_mat, T, n, m);
    del_k = int_matrix(A_mat,T)*arm_dynamics(X(:,k),U(:,k-1)) + X(:,k) -Phi_mat*X(:,k)-Gamma_mat*U(:,k-1);
    Del_k = repmat(del_k,p,1);
    Zs_k = repmat(Xs(:,k), p, 1);
    Us_k = repmat(Us(:,k), p, 1);
    [Sx,Su,Sdel,H,F] = MPC_matrices(Phi_mat,Gamma_mat,p,Wx,Wu,X(:,k),Del_k,Zs_k, Us_k);
    A_ineq = A_inequality(Su);
    b_ineq = b_inequality(X(:,k),U(:,k-1),Del_k,Sx,Sdel);
    Ufk = quadprog(H,F,A_ineq,b_ineq,Aeq,beq,[],[], Ufk_prev,options);
    U(:,k) = Lambda_0'*Ufk;
    [~,Y] = ode45(@(t,Z) arm_dynamics(Z,U(:,k)),[0,T],X(:,k),options_ode); 
    X(:,k+1) = (Y(end,:))';
end

%% Plots
figure;
subplot(2, 2, 1); % Subplot for q1
plot(tspan, X(1, :), 'b', 'LineWidth', 1.5); % Plot q1
hold on;
plot(tspan, Xs(1,:), 'r--', 'LineWidth', 1.5); % Plot setpoint for q1
xlabel('Time (s)');
ylabel('q1 (rad)');
title('State q1 and Setpoint');
legend('q1', 'Setpoint');
grid on;

subplot(2, 2, 2); % Subplot for q2
plot(tspan, X(2, :), 'b', 'LineWidth', 1.5); % Plot q2
hold on;
plot(tspan, Xs(2,:), 'r--', 'LineWidth', 1.5); % Plot setpoint for q2
xlabel('Time (s)');
ylabel('q2 (rad)');
title('State q2 and Setpoint');
legend('q2', 'Setpoint');
grid on;

subplot(2, 2, 3); % Subplot for u1
stairs(tspan(1:end-1), U(1, :), 'g', 'LineWidth', 1.5); % Plot u1
xlabel('Time (s)');
ylabel('u1 (Nm)');
title('Control Input u1');
grid on;

subplot(2, 2, 4); % Subplot for u2
stairs(tspan(1:end-1), U(2, :), 'g', 'LineWidth', 1.5); % Plot u2
xlabel('Time (s)');
ylabel('u2 (Nm)');
title('Control Input u2');
grid on;

% Compute end-effector trajectory
trajectory = [l1*cos(X(1,:)) + l2*cos(X(1,:)+X(2,:)); ...
              l1*sin(X(1,:)) + l2*sin(X(1,:)+X(2,:))];
reference_traj = [l1*cos(Xs(1,:)) + l2*cos(Xs(1,:)+Xs(2,:)); ...
              l1*sin(Xs(1,:)) + l2*sin(Xs(1,:)+Xs(2,:))];

% Animation
figure; axis equal; grid on; hold on;
xlim([-1.1*(l1+l2), 1.1*(l1+l2)]); ylim([-1.1*(l1+l2), 1.1*(l1+l2)]);
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('Manipulator Animation with End-Effector Trajectory');
plot(reference_traj(1,:), reference_traj(2,:));

% Initialize plot handles
h_links = plot([0, 0], [0, 0], 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 2);
h_trajectory = plot(trajectory(:,1), trajectory(:,2), 'b-', 'LineWidth', 1.5);
h_end_effector = plot(0, 0, 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2);

for i = 1:length(tspan)
    q = X(:,i);
    x1 = l1*cos(q(1)); y1 = l1*sin(q(1));
    x2 = x1 + l2*cos(q(1)+q(2)); y2 = y1 + l2*sin(q(1)+q(2));
    
    % Update manipulator links
    set(h_links, 'XData', [0, x1, x2], 'YData', [0, y1, y2]);
    
    % Update end-effector position
    set(h_end_effector, 'XData', x2, 'YData', y2);
    
    % Update trajectory (up to current point)
    set(h_trajectory, 'XData', trajectory(1,1:i), 'YData', trajectory(2,1:i));
    
    pause(0.1);
end

%% Helper Functions

function Gamma = Gamma_matrix(A_matrix, B, T, n, m)
    if det(A_matrix) > 1e-10
        Gamma = (expm(A_matrix * T) - eye(n)) * pinv(A_matrix) * B;
    else
        ode_func = @(tau, x) reshape(expm(A_matrix * tau) * B, [], 1);
        [~, Y] = ode45(ode_func, [0, T], zeros(n * m, 1));
        Gamma = reshape(Y(end, :)', n, m);
    end
end

function int_mat = int_matrix(A_matrix,T)
    n = size(A_matrix,1);
    if det(A_matrix) > 1e-10
        int_mat = (expm(A_matrix * T) - eye(n)) * pinv(A_matrix);
    else
        ode_func = @(tau, x) reshape(expm(A_matrix * tau), [], 1);
        [~, Y] = ode45(ode_func, [0, T], zeros(n * n, 1));
        int_mat = reshape(Y(end, :)', n, n);
    end
end

function Us = find_steady_state_input_trajectory(Xs, m, params)
    % Xs: (nx x N) matrix where each column is a setpoint at a given time.
    % params: struct containing physical parameters (m1, m2, l1, l2, I1, I2, g)
    % Us: (nu x N) matrix where each column is the computed steady-state input.

    [nx, N] = size(Xs);  % N = number of time steps
    Us = zeros(m, N);    % Preallocate for efficiency
    
    % Extract parameters
    m1 = params.m1; m2 = params.m2; 
    l1 = params.l1; l2 = params.l2;
    I1 = params.I1; I2 = params.I2; g = params.g;
    
    % Initial guess for the first step (can be zeros or another guess)
    us_initial_guess = zeros(m, 1);  
    
    % Loop over each setpoint in Xs
    for k = 1:N
        q = Xs(1:2, k);  % Current joint positions
        q_dot = Xs(3:4, k);  % Current joint velocities
        
        % Define residual function for fsolve
        residual = @(us) [
            q_dot(1);  % Steady state requires q_dot = 0
            q_dot(2);
            H_mat(q) \ (us - Cor_mat(q, q_dot)*q_dot - G_mat(q))  % Acceleration should be 0
        ];
        
        % Solve for us at this step
        options = optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');
        us = fsolve(residual, us_initial_guess, options);
        
        % Store the result
        Us(:, k) = us;
        
        % Update initial guess for the next step (warm-start)
        us_initial_guess = us;
    end
end

% Local functions for dynamics calculations
function Hq = H_mat(q)
    % Inertia matrix calculation
    persistent m1 m2 l1 l2 I1 I2 g
    if isempty(m1)
        % Initialize parameters if not set
        m1 = 5; m2 = 4; l1 = 1; l2 = 0.8; 
        I1 = m1*l1^2/3; I2 = m2*l2^2/3; g = 0;
    end
    Hq = [I1 + m1*l1^2/4 + I2 + m2*(l2^2/4 + l1^2 + l1*l2*cos(q(2))), ...
          I2 + m2*(l2^2/4 + l1*l2*cos(q(2))/2); 
         I2 + m2*(l2^2/4 + l1*l2*cos(q(2))/2), ...
         I2 + m2*l2^2/4];
end

function Cq = Cor_mat(q, q_dot)
    % Coriolis matrix calculation
    persistent m2 l1 l2
    if isempty(m2)
        m2 = 4; l1 = 1; l2 = 0.8;
    end
    Cq = [-m2*l1*l2*sin(q(2))*q_dot(2), -m2*l1*l2*sin(q(2))*(q_dot(1)+q_dot(2))/2; 
          m2*l1*l2*sin(q(2))*q_dot(1)/2, 0];
end

function Gq = G_mat(q)
    % Gravity vector calculation
    persistent m1 m2 l1 l2 g
    if isempty(m1)
        m1 = 5; m2 = 4; l1 = 1; l2 = 0.8; g = 0;
    end
    Gq = [m1*g*l1*cos(q(1))/2 + m2*g*(l1*cos(q(1)) + l2/2*cos(q(1)+q(2)));
          m2*g*l2/2*cos(q(1)+q(2))];
end

function [Sx,Su,Sdel,H,F] = MPC_matrices(Phi_mat,Gamma_mat,p,Wx,Wu,x_k,Del_k,Zs,Us)
    n = size(Phi_mat,1); m = size(Gamma_mat,2);
    Sx = zeros(n*p,n); Su = zeros(n*p, m*p); Sdel = zeros(n*p,n*p); % Matrices in equality constraint, X = Sx*x(k) + Su*Uf,k + Sdelta*Deltak
    for i = 1:p
        Sx(n*(i-1)+1:n*i,1:n) = Phi_mat^i; 
        for j = 1:i
            Su(n*(i-1)+1:n*i,m*(i-j)+1:m*(i-j)+m) = Phi_mat^(j-1)*Gamma_mat;
            Sdel(n*(i-1)+1:n*i,n*(i-j)+1:n*(i-j)+n) = Phi_mat^(j-1);
        end
    end
    H = 2*(Su'*Wx*Su + Wu); H = (H + H')/2;
    F = 2*((Wx*Su)'*(Sx*x_k+Sdel*Del_k-Zs)-Wu'*Us);
end