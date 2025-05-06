clear; clc; close all;

%% Physical Parameters for the Quadruple Tank Problem

params.A1 = 28; params.A2 = 32; params.A3 = 28; params.A4 = 32;
params.a1 = 0.071; params.a2 = 0.057; params.a3 = 0.071; params.a4 = 0.057;
params.g = 981;
params.k1 = 3.14; params.k2 = 3.29;
params.gamma1 = 0.43; params.gamma2 = 0.34; params.gamma3 = 0.4;

%% MPC parameters and constraints

T = 4; % Sampling Time (in s)
p = 50; % Prediction Horizon   
q = 20; % Control Horizon
% Steady state Conditions for the above parameters with no disturbance
Us = [3.15;3.15]; Ds = 0; Xs = [12.44;13.17;4.73;4.98]; 
R = Xs; % Setpoint to reach
n = 4; m = 2; % No. of states and control inputs respectively 
wx = eye(n); % State Weights matrix
wu = 0.1*eye(m); % Control weights matrix
U_L = [0;0]; U_H = [5;5]; % Control input constraints
u_L = U_L - Us; u_H = U_H - Us; % For control perturbation
delU_L = [-3;-3]; delU_H = [3;3]; % Control change constraints
X_L = [10;10;2;2]; X_H = [16;16;8;8]; % State constraints
x_L = X_L - Xs; x_H = X_H - Xs; % For state perturbation
load("Continuous_time_linear_model_without_disturbance.mat")
H_mat = [0;0;params.gamma3;(1-params.gamma3)/params.A4]; % Not really needed, just for completeness
Phi = expm(A_mat*T); 
Gamma_u = (Phi-eye(size(Phi)))*pinv(A_mat)*B_mat;
Gamma_d = (Phi-eye(size(Phi)))*pinv(A_mat)*H_mat ; % Not really needed, just for completeness

%% MPC matrices

Wx = kron(eye(p),wx); % using kronecker product we can stack wx in block diagonal matrix p times
Wu = kron(eye(p),wu); % same for wu
Sx = zeros(n*p,n); Su = zeros(n*p, m*p); % Matrices in equality constraint, X = Sx*x(k) + Su*Uf,k
for i = 1:p
    Sx(n*(i-1)+1:n*i,1:n) = Phi^i; 
    for j = 1:i
        Su(n*(i-1)+1:n*i,m*(i-j)+1:m*(i-j)+m) = Phi^(j-1)*Gamma_u;
    end
end
Lambda = zeros(m*p,m*p); Lambda_0 = [eye(m);zeros((p-1)*m,m)]; % Matrices in change of control input inequalities
Lambda(1:m,1:m) = eye(m);
for i = 2:p
    Lambda((i-1)*m + 1:i*m,(i-2)*m + 1:(i-1)*m) = -eye(m);
    Lambda((i-1)*m + 1:i*m,(i-1)*m + 1:i*m) = eye(m);
end
Utilde_L = repmat(u_L,p,1); Utilde_H = repmat(u_H,p,1); % Control constraints stacked
Xtilde_L = repmat(x_L,p,1); Xtilde_H = repmat(x_H,p,1); % State constraints stacked
DeltaU_L = repmat(delU_L,p,1); DeltaU_H = repmat(delU_H,p,1); % Control change constraints stacked

%% QuadProg Matrices for optimization

% X = quadprog(H,f,A,b,Aeq,beq) attempts to solve the quadratic programming 
% problem: min 0.5*x'*H*x + f'*x   subject to:  A*x <= b, Aeq*x = beq
%           x  
H = 2*(Su'*Wx*Su + Wu); H = (H + H')/2; % Although hessian symmetric, ensuring it numerically
f = @(x_k) 2*(Wx*Su)'*Sx*x_k;
A = [eye(m*p);Lambda;Su;-eye(m*p);-Lambda;-Su];
b = @(x_k,u_kminus1) [Utilde_H;  DeltaU_H + Lambda_0*u_kminus1;   Xtilde_H-Sx*x_k;...
                     -Utilde_L;-(DeltaU_L + Lambda_0*u_kminus1);-(Xtilde_L-Sx*x_k)];
Aeq = zeros(m*(p-q),m*p); 
for i = 1:(p-q)
    Aeq((i-1)*m+1:i*m,(q-1)*m+1:q*m) = eye(m);
    Aeq((i-1)*m+1:i*m,(q+i-1)*m+1:(q+i)*m) = -eye(m);
end    
beq = zeros(m*(p-q),1); % These are equality constraints to satisfy constant u 
                        % after control horizon till end of prediction horizon
%% Simulation Variables

Ns = 150; % Number of Simulation steps

% Initialize all matrices
X_case1 = zeros(n,Ns); X_case2 = zeros(n,Ns); % Store states in both cases
x_case1 = zeros(n,Ns); x_case2 = zeros(n,Ns); % Store state perturbations
U_case1 = zeros(m,Ns-1); U_case2 = zeros(m,Ns-1); % Store inputs in both cases
u_case1 = zeros(m,Ns-1); u_case2 = zeros(m,Ns-1); % Store input perturbations

% Set Initial Conditions
X_case1(:,1) = Xs + 2*ones(n,1); X_case2(:,1) = Xs + 2*ones(n,1);
x_case1(:,1) = 2*ones(n,1); x_case2(:,1) = 2*ones(n,1);
U_case1(:,1) = Us; U_case2(:,1) = Us;
u_case1(:,1) = zeros(m,1); u_case(:,2) = zeros(m,1);

%% Simulation
% First Step Simulation with U(0) = Us
x_case1(:,2) = Phi*x_case1(:,1) + Gamma_u*u_case1(:,1) + Gamma_d*0; 
X_case1(:,2) = Xs + x_case1(:,2); 
f_sys = @(t, X) System_Dynamics_210100059(X_case2(:,1), params, U_case2(:,1), 0); 
[~,Y] = ode45(f_sys,[0,T],X_case2(:,1)); 
X_case2(:,2) = (Y(end,:))'; x_case2(:,2) = X_case2(:,2) - Xs;
options = optimoptions('quadprog', 'Display', 'off');
Ufk_case1_prev = zeros(m*p,1); Ufk_case2_prev = zeros(m*p,1); % Setting initial conditions for quadprog 
fprintf('Starting MPC computation for both cases - \n')
for k = 2:Ns-1
    % To display progress
    if mod(k, round(Ns/5)) == 0
        fprintf('%d / %d iterations done\n', k, Ns);
    elseif k == Ns-1 % Print final message at the end
        fprintf('%d / %d iterations done\n', k+1, Ns);
        fprintf("MPC Computation completed")
    end
    Ufk_case1 = quadprog(H,f(x_case1(:,k)),A,b(x_case1(:,k),u_case1(:,k-1)),Aeq,beq,[],[], Ufk_case1_prev,options);
    Ufk_case2 = quadprog(H,f(x_case2(:,k)),A,b(x_case2(:,k),u_case2(:,k-1)),Aeq,beq,[],[], Ufk_case2_prev,options);
    u_case1(:,k) = Lambda_0'*Ufk_case1; U_case1(:,k) = u_case1(:,k) + Us;
    u_case2(:,k) = Lambda_0'*Ufk_case2; U_case2(:,k) = u_case2(:,k) + Us;
    % Case 1 - Plant has same linear model
    x_case1(:,k+1) = Phi*x_case1(:,k) + Gamma_u*u_case1(:,k) + Gamma_d*0; 
    X_case1(:,k+1) = Xs + x_case1(:,k+1); 
    % Case 2 - Plant has non-linear model
    f_sys = @(t, X) System_Dynamics_210100059(X, params, U_case2(:,k), 0); 
    [~,Y] = ode45(f_sys,[0,T],X_case2(:,k)); 
    X_case2(:,k+1) = (Y(end,:))'; x_case2(:,k+1) = X_case2(:,k+1) - Xs;

    Ufk_case1_prev = Ufk_case1; Ufk_case2_prev = Ufk_case2; % Setting initial value for quadprog for next loop
end   
pause(1)
SSE_case1 = sum((X_case1 - R).^2, 2); SSE_case2 = sum((X_case2 - R).^2, 2);
SSMV_case1 = sum((U_case1 - Us).^2, 2); SSMV_case2 = sum((U_case2 - Us).^2, 2);
fprintf('\n\n  Comparison of Metrics for Linear and Nonlinear Plants:\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-25s %-20s %-20s\n', 'Metric', 'Linear Plant', 'Nonlinear Plant');
fprintf('------------------------------------------------------------\n');

for i = 1:length(SSE_case1)
    fprintf('SSE error in height %d         %-20.4f %-20.4f\n', i, SSE_case1(i), SSE_case2(i));
end

for i = 1:length(SSMV_case1)
    fprintf('SSMV error in voltage %d       %-20.4f %-20.4f\n', i, SSMV_case1(i), SSMV_case2(i));
end
pause(1.5)

%% Plotting all graphs
close all;
time = 0:T:(Ns-1)*T;
figure;
set(gcf, 'WindowState', 'maximized'); 
sgtitle("Controlled States", "FontSize", 25, "FontWeight", "bold");
subplot(2,2,1)
hold on
grid on
set(gca,"FontSize",15)
xlabel("Time (in s)","FontSize",20); ylabel("Water level (in cm)","FontSize",20);
title("$h_1$","FontSize",25,"Interpreter","latex")
plot(time, X_case1(1, :), "LineWidth", 3,'Color',[0 0 1], 'DisplayName', '$X_1$ - Linear'); 
plot(time, X_case2(1, :), "LineWidth", 1, 'Marker', 'o', 'MarkerSize', 4, 'DisplayName', '$X_1$ - Non Linear'); 
plot(time, R(1) * ones(size(time)), "LineWidth", 2, 'Color', [1, 1, 0], 'DisplayName', 'Setpoint'); 
plot(time, X_L(1) * ones(size(time)), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Lower Bound');
plot(time, X_H(1) * ones(size(time)), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Upper Bound');
legend("FontSize", 12, "interpreter", "latex")
subplot(2,2,2)
hold on
grid on
set(gca,"FontSize",15)
xlabel("Time (in s)","FontSize",20); ylabel("Water level (in cm)","FontSize",20);
title("$h_2$","FontSize",25,"Interpreter","latex")
plot(time, X_case1(2, :), "LineWidth", 3,'Color',[0 0 1], 'DisplayName', '$X_2 - Linear$');  
plot(time, X_case2(2, :), "LineWidth", 1, 'Marker', 'o', 'MarkerSize', 4, 'DisplayName', '$X_2$ - Non Linear')
plot(time, R(2) * ones(size(time)), "LineWidth", 2, 'Color', [1, 1, 0], 'DisplayName', 'Setpoint'); 
plot(time, X_L(2) * ones(size(time)), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Lower Bound');
plot(time, X_H(2) * ones(size(time)), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Upper Bound');
legend("FontSize", 12, "interpreter", "latex")
subplot(2,2,3)
hold on
grid on
set(gca,"FontSize",15)
xlabel("Time (in s)","FontSize",20); ylabel("Water level (in cm)","FontSize",20);
title("$h_3$","FontSize",25,"Interpreter","latex")
plot(time, X_case1(3, :), "LineWidth", 3,'Color',[0 0 1], 'DisplayName', '$X_3 - Linear$');  
plot(time, X_case2(3, :), "LineWidth", 1, 'Marker', 'o', 'MarkerSize', 4, 'DisplayName', '$X_3$ - Non Linear')
plot(time, R(3) * ones(size(time)), "LineWidth", 2, 'Color', [1, 1, 0], 'DisplayName', 'Setpoint'); 
plot(time, X_L(3) * ones(size(time)), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Lower Bound');
plot(time, X_H(3) * ones(size(time)), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Upper Bound');
legend("FontSize", 12, "interpreter", "latex")
subplot(2,2,4)
hold on
grid on
set(gca,"FontSize",15)
xlabel("Time (in s)","FontSize",20); ylabel("Water level (in cm)","FontSize",20);
title("$h_4$","FontSize",25,"Interpreter","latex")
plot(time, X_case1(4, :), "LineWidth", 3,'Color',[0 0 1], 'DisplayName', '$X_4 - Linear$');
plot(time, X_case2(4, :), "LineWidth", 1, 'Marker', 'o', 'MarkerSize', 4, 'DisplayName', '$X_4$ - Non Linear')
plot(time, R(4) * ones(size(time)), "LineWidth", 2, 'Color', [1, 1, 0], 'DisplayName', 'Setpoint'); 
plot(time, X_L(4) * ones(size(time)), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Lower Bound');
plot(time, X_H(4) * ones(size(time)), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Upper Bound');
legend("FontSize", 12, "interpreter", "latex")

figure;
set(gcf, 'WindowState', 'maximized'); % Maximizes the figure window
sgtitle("Manipulated Inputs", "FontSize", 25, "FontWeight", "bold");
subplot(2,1,1)
hold on
grid on
set(gca,"FontSize",15)
xlabel("Time (in s)","FontSize",20); ylabel("Voltage level (in Volts)","FontSize",20);
title("$v_1$","FontSize",25,"Interpreter","latex")
stairs(time(1:end-1),U_case1(1,:),"LineWidth",3,'Color',[0 0 1])
stairs(time(1:end-1),U_case2(1,:),"LineWidth",1, 'Marker', 'o', 'MarkerSize', 4)
plot(time(1:end-1),U_L(1)*ones(length(time)-1,1),"LineWidth",3,"LineStyle","-.")
plot(time(1:end-1),U_H(1)*ones(length(time)-1,1),"LineWidth",3,"LineStyle","-.")
legend("$v_1$ - Linear","$v_1$ - Non-Linear","Lower Bound","Upper Bound","FontSize",12,"interpreter","latex")
subplot(2,1,2)
hold on
grid on
set(gca,"FontSize",15)
xlabel("Time (in s)","FontSize",20); ylabel("Voltage level (in Volts)","FontSize",20);
title("$v_2$","FontSize",25,"Interpreter","latex")
stairs(time(1:end-1),U_case1(2,:),"LineWidth",3,'Color',[0 0 1])
stairs(time(1:end-1),U_case2(2,:),"LineWidth",1, 'Marker', 'o', 'MarkerSize', 4)
plot(time(1:end-1),U_L(2)*ones(length(time)-1,1),"LineWidth",3,"LineStyle","-.")
plot(time(1:end-1),U_H(2)*ones(length(time)-1,1),"LineWidth",3,"LineStyle","-.")
legend("$v_2$ - Linear","$v_2$ - Non-Linear","Lower Bound","Upper Bound","FontSize",12,"interpreter","latex")
