function [X,U,SSE,SSMV] = MPC()
    
    %clear; clc; close all;
    
    %% Physical Parameters for the Quadruple Tank Problem
    
    params.A1 = 28; params.A2 = 32; params.A3 = 28; params.A4 = 32;
    params.a1 = 0.071; params.a2 = 0.057; params.a3 = 0.071; params.a4 = 0.057;
    params.g = 981;
    params.k1 = 3.33; params.k2 = 3.35;
    params.gamma1 = 0.7; params.gamma2 = 0.6; params.gamma3 = 0.4;
    
    %% MPC parameters and constraints
    
    T = 4; % Sampling Time (in s)
    p = 50; % Prediction Horizon   
    q = 20; % Control Horizon
    % Steady state Conditions for the above parameters with no disturbance
    Us = [3;3]; Ds = 0;  Xs = [12.263;12.7831;1.6339;1.409];
    R = Xs; % Setpoint to reach
    n = 4; m = 2; % No. of states and control inputs respectively 
    wx = eye(n); % State Weights matrix
    wu = 0.1*eye(m); % Control weights matrix
    U_L = [0;0]; U_H = [5;5]; % Control input constraints
    u_L = U_L - Us; u_H = U_H - Us; % For control perturbation
    delU_L = [-3;-3]; delU_H = [3;3]; % Control change constraints
    X_L = [10;10;1;1]; X_H = [16;16;8;8]; % State constraints
    x_L = X_L - Xs; x_H = X_H - Xs; % For state perturbation
    load("Continuous_time_linear_model_without_dist.mat")
    Phi = expm(A_mat*T); 
    Gamma_u = (Phi-eye(size(Phi)))*pinv(A_mat)*B_mat;
    %Gamma_d = (Phi-eye(size(Phi)))*pinv(A_mat)*D_mat ; % Not really needed, just for completeness
    
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
    
    Ns = 200; % Number of Simulation steps
    
    % Initialize all matrices
    X = zeros(n,Ns); x = zeros(n,Ns);
    U = zeros(m,Ns-1); u = zeros(m,Ns-1); 
    
    % Set Initial Conditions
    X(:,1) = Xs + 2*ones(n,1); x(:,1) = 2*ones(n,1); 
    U(:,1) = Us; u(:,1) = zeros(m,1);
    
    %% Simulation
    % First Step Simulation with U(0) = Us
    f_sys = @(t, X) System_Dynamics_210100059(X(:,1), params, U(:,1), 0); 
    [~,Y] = ode45(f_sys,[0,T],X(:,1)); 
    X(:,2) = (Y(end,:))'; x(:,2) = X(:,2) - Xs;
    options = optimoptions('quadprog','Algorithm','active-set', 'Display', 'off');
    Ufk_prev = zeros(m*p,1); % Setting initial conditions for quadprog 
    %fprintf('Starting MPC computation - \n')
    for k = 2:Ns-1
        % % To display progress
        % if mod(k, round(Ns/5)) == 0
        %     fprintf('%d / %d iterations done\n', k, Ns);
        % elseif k == Ns-1 % Print final message at the end
        %     fprintf('%d / %d iterations done\n', k+1, Ns);
        %     fprintf("MPC Computation completed")
        % end
        Ufk = quadprog(H,f(x(:,k)),A,b(x(:,k),u(:,k-1)),Aeq,beq,[],[], Ufk_prev,options);
        u(:,k) = Lambda_0'*Ufk; U(:,k) = u(:,k) + Us;
        f_sys = @(t, X) System_Dynamics_210100059(X, params, U(:,k), 0); 
        [~,Y] = ode45(f_sys,[0,T],X(:,k)); 
        X(:,k+1) = (Y(end,:))'; x(:,k+1) = X(:,k+1) - Xs;
        Ufk_prev = Ufk; % Setting initial value for quadprog for next loop
    end   
    
    SSE = sum((X - R).^2, 2);
    SSMV = sum((U - Us).^2, 2); 

end
% fprintf('\n\n  Metrics for the MPC Controller:\n');
% fprintf('-------------------------------------------\n');
% fprintf('%-25s %-20s\n', 'Metric', 'Value');
% fprintf('-------------------------------------------\n');
% 
% for i = 1:length(SSE)
%     fprintf('SSE error in height %d         %-20.4f\n', i, SSE(i));
% end
% 
% for i = 1:length(SSMV)
%     fprintf('SSMV error in voltage %d       %-20.4f\n', i, SSMV(i));
% end
% 
% %% Plotting
% 
% close all;
% time = 0:T:(Ns-1)*T;
% 
% % Plot Controlled States
% figure;
% set(gcf, 'WindowState', 'maximized'); 
% sgtitle("MPC Controller - Controlled States", "FontSize", 25, "FontWeight", "bold");
% 
% for i = 1:4
%     subplot(2,2,i)
%     hold on
%     grid on
%     set(gca,"FontSize",15)
%     xlabel("Time (in s)","FontSize",20); 
%     ylabel("Water level (in cm)","FontSize",20);
%     title(sprintf("$h_%d$", i), "FontSize", 25, "Interpreter", "latex")
% 
%     plot(time, X(i, :), "LineWidth", 3, 'DisplayName', 'PPC');  
%     plot(time, R(i) * ones(size(time)), "LineWidth", 2, 'DisplayName', 'Setpoint'); 
%     plot(time, X_L(i) * ones(size(time)), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Lower Bound');
%     plot(time, X_H(i) * ones(size(time)), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Upper Bound');
% 
%     legend("FontSize", 12, "interpreter", "latex")
% end
% 
% % Plot Manipulated Inputs
% figure;
% set(gcf, 'WindowState', 'maximized'); 
% sgtitle("MPC Controller - Manipulated Inputs", "FontSize", 25, "FontWeight", "bold");
% 
% for i = 1:2
%     subplot(2,1,i)
%     hold on
%     grid on
%     set(gca,"FontSize",15)
%     xlabel("Time (in s)","FontSize",20); 
%     ylabel("Voltage level (in Volts)","FontSize",20);
%     title(sprintf("$v_%d$", i), "FontSize", 25, "Interpreter", "latex")
% 
%     stairs(time(1:end-1), U(i,:), "LineWidth", 3, 'Color', [0 0 1], 'DisplayName', 'PPC');
%     plot(time(1:end-1), U_L(i) * ones(length(time)-1,1), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Lower Bound');
%     plot(time(1:end-1), U_H(i) * ones(length(time)-1,1), "LineWidth", 3, "LineStyle", "-.", 'DisplayName', 'Upper Bound');
% 
%     legend("FontSize", 12, "interpreter", "latex")
% end
