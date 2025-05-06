clc; clear; close all;

%% System Parameters

m1 = 1; % Mass of link 1
m2 = 1; % Mass of link 2
m3 = 1; % Mass of link 3
l1 = 1; % Length of link 1
l2 = 1; % Length of link 2
l3 = 1; % Length of link 3
g = 0; % Acceleration due to gravity 
I1 = (m1*l1^2)/12; % Moment of inertia of link 1 about COM 
I2 = (m2*l2^2)/12; % Moment of inertia of link 2 about COM
I3 = (m3*l3^2)/12; % Moment of inertia of link 3 about COM
lc1 = (l1/2); % Location of COM of link 1
lc2 = (l2/2); % Location of COM of link 2
lc3 = (l3/2); % Location of COM of link 3
params = struct('m1', m1,'m2', m2, 'm3', m3,'l1', l1, 'l2', l2, 'l3', l3, 'g', g, ...
    'I1', I1, 'I2', I2, 'I3', I3, 'lc1', lc1, 'lc2', lc2, 'lc3', lc3 ...
     );

T = 0.04; % Time step
Ns = 1000; % No. of Time Steps
Tmax = Ns*T; % Maximum time

tspan = 0:T:Tmax; % Time duration for the simulation

reach_vec = zeros(Ns-1,1); % Vector storing rank of reachability matrix at each time step
obs_vec = zeros(Ns-1,1); % Vector storing rank of observability matrix at each time step
stab_vec = zeros(Ns-1,1); % Vector storing spectral radius of Phi at each time step

%% MPC Parameters

Np = 30; % Prediction Horizon   
q = 20; % Control Horizon
n = 8; % No. of states (including base and manipulator states)
m = 4; % No. of control inputs 
p = 8; % No. of outputs (x and y position of end effector)

x0 = zeros(n,1); % Initial state

% State and Control Weights
wx = eye(n); % State Weights matrix
wu = 0.01 * eye(m); % Control weights matrix

% Observer weights
Q = eye(n); R_mat = 10*eye(p);

% Control input constraints
U_L = -5 * ones(m, 1); 
U_H = 5 * ones(m, 1); 

% Control change constraints
delU_L = -0.5 * ones(m, 1); 
delU_H = 0.5 * ones(m, 1); 

% State constraints
X_L = -4*pi * ones(n, 1); 
X_H = 4*pi * ones(n, 1); 

% MPC Constant matrices
Wx = kron(eye(Np), wx); % State weights matrix
Wu = kron(eye(Np), wu); % Control weights matrix
Lambda = zeros(m * Np, m * Np); 
Lambda_0 = [eye(m); zeros((Np - 1) * m, m)]; % Matrices in change of control input inequalities
Lambda(1:m, 1:m) = eye(m);
for k = 2:Np
    Lambda((k - 1) * m + 1:k * m, (k - 2) * m + 1:(k - 1) * m) = -eye(m);
    Lambda((k - 1) * m + 1:k * m, (k - 1) * m + 1:k * m) = eye(m);
end
Utilde_L = repmat(U_L, Np, 1); 
Utilde_H = repmat(U_H, Np, 1); 
Xtilde_L = repmat(X_L, Np, 1); 
Xtilde_H = repmat(X_H, Np, 1); 
DeltaU_L = repmat(delU_L, Np, 1); 
DeltaU_H = repmat(delU_H, Np, 1); 

Aeq = zeros(m * (Np - q), m * Np); 
for i = 1:(Np - q)
    Aeq((i - 1) * m + 1:i * m, (q - 1) * m + 1:q * m) = eye(m);
    Aeq((i - 1) * m + 1:i * m, (q + i - 1) * m + 1:(q + i) * m) = -eye(m);
end    
beq = zeros(m * (Np - q), 1);

A_inequality = @(Su) [eye(m * Np); Lambda; Su; -eye(m * Np); -Lambda; -Su];
b_inequality = @(x_k, u_kminus1, Del_k, Sx, Sdel) [Utilde_H; DeltaU_H + Lambda_0 * u_kminus1; Xtilde_H - Sx * x_k - Sdel * Del_k; ...
                                         -Utilde_L; -DeltaU_L - Lambda_0 * u_kminus1; -Xtilde_L + Sx * x_k + Sdel * Del_k];

%% MPC simulation
X = zeros(n, Ns + 1); X(:, 1) = x0;X_pred = zeros(n,Ns + 1); X_est = zeros(n,Ns + 1);
U = zeros(m, Ns); Y = zeros(p, Ns + 1);
Xs = [pi/3*ones(size(tspan));-pi/3*ones(size(tspan));-5*ones(size(tspan));0.5*ones(size(tspan));0*ones(size(tspan));0*ones(size(tspan));0*ones(size(tspan));0*ones(size(tspan))];
%Xs = [pi*cos(0.2*tspan);pi*sin(0.2*tspan);-6*cos(0.2*tspan);4*sin(0.2*tspan);...
%    -pi*0.2*sin(0.2*tspan);pi*0.2*cos(0.2*tspan);6*0.2*sin(0.2*tspan);4*0.2*cos(0.2*tspan)];
Us = find_steady_state_input_trajectory(Xs, m, params);
X(:, 1) = Xs(:,1);
X_pred(:, 1) = X(:,1)+0.1*ones(n,1); X_est(:,1) = X(:,1)+0.1*ones(n,1);
Ufk_prev = zeros(m * Np, 1);
options_ode = odeset('abstol', 1e-9, 'reltol', 1e-9);
options = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off');
% First step simulation with U(0) = 0
A_mat = A_matrix(X_est(:, 1), U(:, 1),params);
B_mat = B_matrix(X_est(:, 1), U(:, 1),params);
%C_mat = C_matrix(X_est(:, 1), params);
C_mat = eye(8);
Phi_mat = expm(A_mat * T);
A_const = A_matrix(Xs(:, 1), Us(:, 1),params);
B_const = B_matrix(Xs(:, 1), Us(:, 1),params);
Phi_const = expm(T*A_const);
Gamma_const = Gamma_matrix(A_const, B_const, T, n, m);
Gamma_mat = Gamma_matrix(A_mat, B_mat, T, n, m);
p_des = [0.9;0.6;0.8;0.85;0.9;0.95;0.6;0.5];
%L_const = observer(Phi_const,C_mat,Q,R_mat);
L_const = place(Phi_const',C_mat',p_des);
[~,val] = ode45(@(t,Z) dyn(Z,U(:,1),params),[0,T],X(:,1),options_ode); 
X(:,2) = (val(end,:))';
% Measurement
Y(:,1:2) = end_effec(X(:,1:2), params);
% Prediction step :
X_pred(:,2) = Xs(:,1)+Phi_const*(X_est(:,1)-Xs(:,1))+ Gamma_const*(U(:,1)-Us(:,1));%dyn(X_est(:,1),U(:,k),params);
% Update Step
X_est(:,2) = X_pred(:,1) + L_const*(Y(:,2)-C_mat*(X_pred(:,2)-Xs(:,2)));


fprintf("Starting MPC Computation : \n")
for k = 2:Ns
    % To display progress
    if mod(k, round(Ns / 5)) == 0
        fprintf('%d / %d iterations done\n', k, Ns);
    end
    
    % Matrix calculations using estimated states and control
    A_mat = A_matrix(X_est(:, k), U(:, k-1),params);
    B_mat = B_matrix(X_est(:, k), U(:, k-1),params);
    %C_mat = C_matrix(X_est(:,k), params);
    C_mat = eye(8);

    Phi_mat = expm(A_mat * T);
    Gamma_mat = Gamma_matrix(A_mat, B_mat, T, n, m);
    int_mat = int_matrix(A_mat,T);
    %L = observer(Phi_mat,C_mat,Q,R_mat);

    % reach_vec(k) = rank(ctrb(Phi_mat,Gamma_mat));
    % obs_vec(k) = rank(obsv(Phi_mat,C_mat));
    % stab_vec(k) = max(abs(eig(Phi_mat)));

    % % MPC Control Calculation

    del_k = int_mat* dyn(X_est(:,k),U(:,k-1),params)+ X_est(:,k) -Phi_mat*X_est(:,k)-Gamma_mat*U(:,k-1);
    Del_k = repmat(del_k, Np, 1);
    Zs_k = repmat(Xs(:,k), Np, 1);
    Us_k = repmat(Us(:,k), Np, 1);

    [Sx, Su, Sdel, H, F] = MPC_matrices(Phi_mat, Gamma_mat, Np, Wx, Wu, X_est(:, k), Del_k, Zs_k, Us_k);
    A_ineq = A_inequality(Su);
    b_ineq = b_inequality(X_est(:, k), U(:, k-1), Del_k, Sx, Sdel);

    Ufk = quadprog(H, F, A_ineq, b_ineq, Aeq, beq, [], [], Ufk_prev, options);
    U(:, k) = Lambda_0' * Ufk;
    %U(:, k) = [sin(k);cos(k);sin(0.2*k);cos(0.2*k)];

    % Simulate the system forward
    [~, val] =  ode45(@(t,Z) dyn(Z,U(:,k),params),[0,T],X(:,k),options_ode);
    X(:, k + 1) = (val(end, :))';

    % Measurement
    Y(:,k+1) = end_effec(X(:,k+1), params);
    % Prediction step :
    X_pred(:,k+1) = Xs(:,k)+Phi_const*(X_est(:,k)-Xs(:,k))+ Gamma_const*(U(:,k)-Us(:,k));%dyn(X_est(:,k),U(:,k),params);
    % Update Step
    X_est(:,k+1) = X_pred(:,k+1) + L_const*(Y(:,k+1)-C_mat*(X_pred(:,k+1)-Xs(:,k)));
    
    % Update previous control input
    Ufk_prev = Ufk;
end

% %% Animation
% figure;
% set(gcf,"WindowState","maximized")
% axis equal; grid on;
% hold on;
% xlim([-10, 10]); ylim([-10, 10]);
% xlabel('X Position (m)'); ylabel('Y Position (m)');
% title('Loco-Manipulator Animation');
% %h_end_effector = plot(0, 0, 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2);
% %traj_ef = 
% 
% for i = 1:5:length(tspan)
%     x1 = X(3, i) + l1 * cos(X(1, i));
%     y1 = X(4, i) + l1 * sin(X(1, i));
%     x2 = x1 + l2 * cos(X(1, i) + X(2, i));
%     y2 = y1 + l2 * sin(X(1, i) + X(2, i));
% 
%     %scatter(X(3, i), X(4, i), 1000, "black", "filled")
%     %hold on
%     plot([X(3, i), x1], [X(4, i), y1], 'r', 'LineWidth', 10)
%     hold on
%     plot([x1, x2], [y1, y2], 'b-o', 'LineWidth', 2)
%     plot(Xs(3,:),Xs(4,:),"LineStyle","--","LineWidth",2);
% 
% 
%     grid on, set(gca, 'FontSize', 15)
%     axis([-10, 10, -10, 10])
%     axis equal
%     xlabel('$x$', 'Interpreter', 'latex')
%     ylabel('$y$', 'Interpreter', 'latex')
%     title('Animation', 'Interpreter', 'latex')
%     time = Tmax * (i - 1) / (length(tspan) - 1);
%     subtitle("Time = " + time + " sec", 'Interpreter', 'latex')
%     hold off
%     %set(h_trajectory, 'XData', trajectory(1,1:i), 'YData', trajectory(2,1:i));
% 
%     pause(0.01)
% end

%% Plots
figure;
subplot(2, 2, 1); % Subplot for q1
plot(tspan, X(1, 1:end), 'b', 'LineWidth', 1.5); % Plot q1
hold on;
plot(tspan,Xs(1,:),'r','LineWidth',1.5,"LineStyle","--")
xlabel('Time (s)');
ylabel('q1 (rad)');
title('State q1');
grid on;

subplot(2, 2, 2); % Subplot for q2
plot(tspan, X(2, 1:end), 'b', 'LineWidth', 1.5); % Plot q2
hold on;
plot(tspan,Xs(2,:),'r','LineWidth',1.5,"LineStyle","--")
xlabel('Time (s)');
ylabel('q2 (rad)');
title('State q2');
grid on;

subplot(2, 2, 3); % Subplot for u1
stairs(tspan(1:end-1), U(1, :), 'g', 'LineWidth', 1.5); % Plot u1
xlabel('Time (s)');
ylabel('\tau_1 (Nm)');
title('Control Input \tau_1');
grid on;

subplot(2, 2, 4); % Subplot for u2
stairs(tspan(1:end-1), U(2, :), 'g', 'LineWidth', 1.5); % Plot u2
xlabel('Time (s)');
ylabel('\tau_2 (Nm)');
title('Control Input \tau_2');
grid on;

figure;
subplot(2, 2, 1); % Subplot for x
plot(tspan, X(3, 1:end), 'b', 'LineWidth', 1.5); % Plot x
hold on;
plot(tspan,Xs(3,:),'r','LineWidth',1.5,"LineStyle","--")
xlabel('Time (s)');
ylabel('x (m)');
title('State x');
grid on;

subplot(2, 2, 2); % Subplot for y
plot(tspan, X(4, 1:end), 'b', 'LineWidth', 1.5); % Plot y
hold on;
plot(tspan,Xs(4,:),'r','LineWidth',1.5,"LineStyle","--")
xlabel('Time (s)');
ylabel('y m');
title('State y');
grid on;

subplot(2, 2, 3); % Subplot for u3
stairs(tspan(1:end-1), U(3, :), 'g', 'LineWidth', 1.5); % Plot u3
xlabel('Time (s)');
ylabel('Fx (N)');
title('Control Input Fx');
grid on;

subplot(2, 2, 4); % Subplot for u4
stairs(tspan(1:end-1), U(4, :), 'g', 'LineWidth', 1.5); % Plot u4
xlabel('Time (s)');
ylabel('Fy (Nm)');
title('Control Input Fy');
grid on;

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
    % params: Additional parameters for dyn(x,u,params).
    % Us: (nu x N) matrix where each column is the computed steady-state input.

    [nx, N] = size(Xs);  % N = number of time steps
    Us = zeros(m, N);    % Preallocate for efficiency
    
    % Initial guess for the first step (can be zeros or another guess)
    us_initial_guess = zeros(m, 1);  
    
    % Loop over each setpoint in Xs
    for k = 1:N
        xs = Xs(:, k);  % Current setpoint
        
        % Define residual function for fsolve
        residual = @(us) dyn(xs, us, params);
        
        % Solve for us at this step
        options = optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');  % Suppress output
        us = fsolve(residual, us_initial_guess, options);
        
        % Store the result
        Us(:, k) = us;
        
        % Update initial guess for the next step (warm-start)
        us_initial_guess = us;
    end
end

function L = observer(Phi,C_mat,Q,R)
    [P_steady,~, ~] = idare(Phi', C_mat', Q, R, [], []);
    L = P_steady * C_mat' * pinv(C_mat * P_steady * C_mat' + R);
end

function Y = end_effec(X, params)
    % l1 = params.l1; l2 = params.l2;
    % th1 = X(1,:); th2 = X(2,:); px = X(3,:); py = X(4,:);
    % Y(1,:) = px + l1*cos(th1) + l2*cos(th1+th2);
    % Y(2,:) = py + l1*sin(th1) + l2*sin(th1+th2);
    Y = eye(8)*X;
end
