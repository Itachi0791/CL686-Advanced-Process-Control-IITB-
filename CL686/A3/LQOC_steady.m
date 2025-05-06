function [X,U,SSE,SSMV] = LQOC_steady()
    %clc; clear; close all;
    
    %% Physical Parameters for the Quadruple Tank Problem
    
    params.A1 = 28; params.A2 = 32; params.A3 = 28; params.A4 = 32;
    params.a1 = 0.071; params.a2 = 0.057; params.a3 = 0.071; params.a4 = 0.057;
    params.g = 981;
    params.k1 = 3.33; params.k2 = 3.35;
    params.gamma1 = 0.7; params.gamma2 = 0.6; params.gamma3 = 0.4;
    n = 4; m=2;
    
    %% Simulation Conditions
    
    Us = [3;3]; Ds = 0; Xs = [12.263;12.7831;1.6339;1.409];
    T = 4;
    Ns = 200;
    X0 = Xs + [2;2;2;2];
    R = Xs;
    X_L = [10;10;1;1]; X_H = [16;16;8;8];
    U_L = [0;0]; U_H = [5;5]; 
    delU_L = [-3;-3]; delU_H = [3;3];
    load("Continuous_time_linear_model_without_dist.mat")
    Phi = expm(A_mat*T); 
    Gamma_u = (Phi-eye(size(Phi)))*pinv(A_mat)*B_mat;
    
    %% Controller Design
    Wx = eye(n); Wu = 0.1*eye(m);
    [G_inf,~,~] = dlqr(Phi,Gamma_u,Wx,Wu);
    
    %% Simulation
    
    X = zeros(n,Ns); X(:,1) = X0; x = zeros(n,Ns); x(:,1) = X0-Xs;
    U = zeros(m,Ns-1); U_prev = Us;
    u = zeros(m,Ns-1); u_prev = 0; % Assuming U = Us before simulation starts
    for i=1:Ns-1
        u_i = -G_inf*x(:,i); U_i = u_i + Us;
        delU_i = max(min(U_i-U_prev,delU_H),delU_L); % Sets deltaU between delUH and delUL
        U(:,i) = max(min(U_prev + delU_i, U_H), U_L); % Sets between U_H and U_L 
        f_sys = @(t, X) System_Dynamics_210100059(X, params, U(:,i), 0); 
        [~,Y] = ode45(f_sys,[0,T],X(:,i)); 
        X(:,i+1) = (Y(end,:))';
        u_prev = u_i; x(:,i+1) = X(:,i+1)-Xs;
    end
   % pause(0.5)
    SSE = sum((X - R).^2, 2);
    SSMV = sum((U - Us).^2, 2); 
    % 
    % fprintf('\n\n  Metrics for the LQOC (Steady):\n');
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
    % pause(1)
    % %% Plotting
    % 
    % close all;
    % time = 0:T:(Ns-1)*T;
    % 
    % % Plot Controlled States
    % figure;
    % set(gcf, 'WindowState', 'maximized'); 
    % sgtitle("LQOC (Steady) - Controlled States", "FontSize", 25, "FontWeight", "bold");
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
    % sgtitle("LQOC (Steady) - Manipulated Inputs", "FontSize", 25, "FontWeight", "bold");
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
end