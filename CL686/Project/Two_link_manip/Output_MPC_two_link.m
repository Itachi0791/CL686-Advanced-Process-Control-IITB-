clc; clear; close all;

% Physical parameters
m1 = 5; m2 = 4; l1 = 1; l2 = 0.8; I1 = m1*l1^2/3; I2 = m2*l2^2/3; g = 0;

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
Np = 80; % Prediction Horizon   
Nq = 50; % Control Horizon
n = 4; m = 2; p =2; % No. of states, control inputs and outputs respectively 
%Xs = [3*pi/4*sin(0.2*tspan);pi/8-pi/3*cos(0.2*tspan);3*pi/4*0.2*cos(0.2*tspan);pi/3*0.2*sin(0.2*tspan)];  
Ys = [1.5*ones(size(tspan));0.1*ones(size(tspan))]; % Output Setpoint
%Us = G(Xs(1:2)); % Control Setpoint
wy = eye(p); %  Weights matrix
wu = 0.001*eye(m); % Control weights matrix
U_L = -20*ones(m,1); U_H = 20*ones(m,1); % Control input constraints
delU_L = -5*ones(m,1); delU_H = 5*ones(m,1); % Control change constraints
X_L = [-pi;-pi;-5;-5]; X_H = [pi;pi;5;5]; % State constraints
Y_L = [-(l1+l2);-(l1+l2)]; Y_H = [l1+l2;l1+l2]; % Output Constraints

%% MPC Constant matrices

Wy = kron(eye(Np),wy); % using kronecker product we can stack wx in block diagonal matrix p times
Wu = kron(eye(Np),wu); % same for wu
Lambda = zeros(m*Np,m*Np); Lambda_0 = [eye(m);zeros((Np-1)*m,m)]; % Matrices in change of control input inequalities
Lambda(1:m,1:m) = eye(m);
for k = 2:Np
    Lambda((k-1)*m + 1:k*m,(k-2)*m + 1:(k-1)*m) = -eye(m);
    Lambda((k-1)*m + 1:k*m,(k-1)*m + 1:k*m) = eye(m);
end
Utilde_L = repmat(U_L,Np,1); Utilde_H = repmat(U_H,Np,1); % Control constraints stacked
Xtilde_L = repmat(X_L,Np,1); Xtilde_H = repmat(X_H,Np,1); % State constraints stacked
Ytilde_L = repmat(Y_L,Np,1); Ytilde_H = repmat(Y_H,Np,1); % Output constraints stacked
DeltaU_L = repmat(delU_L,Np,1); DeltaU_H = repmat(delU_H,Np,1); % Control change constraints stacked
Aeq = zeros(m*(Np-Nq),m*Np); 
for i = 1:(Np-Nq)
    Aeq((i-1)*m+1:i*m,(Nq-1)*m+1:Nq*m) = eye(m);
    Aeq((i-1)*m+1:i*m,(Nq+i-1)*m+1:(Nq+i)*m) = -eye(m);
end    
beq = zeros(m*(Np-Nq),1);
A_inequality = @(C_big,Su)[eye(m*Np);Lambda;Su;C_big*Su;-eye(m*Np);-Lambda;-Su;-C_big*Su];
b_inequality = @(C_big,x_k,u_kminus1,Del_k,Sx,Sdel) [Utilde_H;DeltaU_H+Lambda_0*u_kminus1;Xtilde_H-Sx*x_k-Sdel*Del_k;Ytilde_H-C_big*(Sx*x_k+Sdel*Del_k);...
                                         -Utilde_L;-DeltaU_L-Lambda_0*u_kminus1;-Xtilde_L+Sx*x_k+Sdel*Del_k;-Ytilde_L+C_big*(Sx*x_k+Sdel*Del_k)];

%% MPC simulation

X = zeros(n,Ns+1); U = zeros(m,length(tspan)-1); 
X0 = [0;0;0;0]; X(:,1) = X0;% U(:,1) = Us;
options_ode = odeset('abstol', 1e-9, 'reltol', 1e-9);
% First step simulation with U(0) = 0
[~,val] = ode45(@(t,Z) arm_dynamics(Z,U(:,1)),[0,T],X(:,1),options_ode); 
X(:,2) = (val(end,:))';
Ufk_prev = zeros(m*Np,1);
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
    C_mat = C_matrix(X(:,k),U(:,k-1));
    Phi_mat = expm(A_mat*T);
    Gamma_mat = Gamma_matrix(A_mat, B_mat, T, n, m);
    del_k = T*arm_dynamics(X(:,k),U(:,k-1)) + X(:,k) -Phi_mat*X(:,k)-Gamma_mat*U(:,k-1);
    Del_k = repmat(del_k,Np,1);
    C_big = kron(eye(Np),C_mat);
    Zs = repmat(Ys(:,k), Np, 1);
    [Sx,Su,Sdel,H,F] = MPC_matrices(Phi_mat,Gamma_mat,C_big,Np,Wy,Wu,X(:,k),Del_k,Zs);
    A_ineq = A_inequality(C_big,Su);
    b_ineq = b_inequality(C_big,X(:,k),U(:,k-1),Del_k,Sx,Sdel);
    Ufk = quadprog(H,F,A_ineq,b_ineq,Aeq,beq,[],[], Ufk_prev,options);
    U(:,k) = Lambda_0'*Ufk;
    [~,val] = ode45(@(t,Z) arm_dynamics(Z,U(:,k)),[0,T],X(:,k),options_ode); 
    X(:,k+1) = (val(end,:))';
end

Y = [l1*cos(X(1,:)) + l2*cos(X(1,:)+X(2,:)); ...
              l1*sin(X(1,:)) + l2*sin(X(1,:)+X(2,:))];

%% Plots
figure;
subplot(2, 2, 1); % Subplot for x
plot(tspan, Y(1, :), 'b', 'LineWidth', 1.5); % Plot x
hold on;
plot(tspan, Ys(1,:), 'r--', 'LineWidth', 1.5); % Plot setpoint for x
xlabel('Time (s)');
ylabel('x (m)');
title('X coordinate');
legend('x', 'Setpoint');
grid on;

subplot(2, 2, 2); % Subplot for q2
plot(tspan, Y(2, :), 'b', 'LineWidth', 1.5); % Plot y
hold on;
plot(tspan, Ys(2,:), 'r--', 'LineWidth', 1.5); % Plot setpoint for y
xlabel('Time (s)');
ylabel('y (m)');
title('Y coordinate');
legend('y', 'Setpoint');
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

% Animation
figure; axis equal; grid on; hold on;
xlim([-1.1*(l1+l2), 1.1*(l1+l2)]); ylim([-1.1*(l1+l2), 1.1*(l1+l2)]);
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('Manipulator Animation with End-Effector Trajectory');
plot(Ys(1,:), Ys(2,:));

% Initialize plot handles
h_links = plot([0, 0], [0, 0], 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 2);
h_trajectory = plot(Y(:,1), Y(:,2), 'b-', 'LineWidth', 1.5);
h_end_effector = plot(0, 0, 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2);

for i = 1:length(tspan)
    q_val = X(:,i);
    x1 = l1*cos(q_val(1)); y1 = l1*sin(q_val(1));
    x2 = x1 + l2*cos(q_val(1)+q_val(2)); y2 = y1 + l2*sin(q_val(1)+q_val(2));
    
    % Update manipulator links
    set(h_links, 'XData', [0, x1, x2], 'YData', [0, y1, y2]);
    
    % Update end-effector position
    set(h_end_effector, 'XData', x2, 'YData', y2);
    
    % Update trajectory (up to current point)
    set(h_trajectory, 'XData', Y(1,1:i), 'YData', Y(2,1:i));
    
    pause(0.05);
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

function [Sx,Su,Sdel,H,F] = MPC_matrices(Phi_mat,Gamma_mat,C_big,p,Wy,Wu,x_k,Del_k,Zs)
    n = size(Phi_mat,1); m = size(Gamma_mat,2);
    Sx = zeros(n*p,n); Su = zeros(n*p, m*p); Sdel = zeros(n*p,n*p); % Matrices in equality constraint, X = Sx*x(k) + Su*Uf,k + Sdelta*Deltak
    for i = 1:p
        Sx(n*(i-1)+1:n*i,1:n) = Phi_mat^i; 
        for j = 1:i
            Su(n*(i-1)+1:n*i,m*(i-j)+1:m*(i-j)+m) = Phi_mat^(j-1)*Gamma_mat;
            Sdel(n*(i-1)+1:n*i,n*(i-j)+1:n*(i-j)+n) = Phi_mat^(j-1);
        end
    end
    H = 2*(Su'*C_big'*Wy*C_big*Su + Wu); H = (H + H')/2;
    F = 2*(Wy*C_big*Su)'*(C_big*(Sx*x_k+Sdel*Del_k)-Zs);
end