clear; clc; close all;
    
%% Physical Parameters for the Quadruple Tank Problem

params.A1 = 28; params.A2 = 32; params.A3 = 28; params.A4 = 32;
params.a1 = 0.071; params.a2 = 0.057; params.a3 = 0.071; params.a4 = 0.057;
params.g = 981;
params.k1 = 3.33; params.k2 = 3.35;
params.gamma1 = 0.7; params.gamma2 = 0.6; params.gamma3 = 0.4;

%% Simulation Parameters
Xs = [12.263;12.7831;1.6339;1.409]; Us = [3;3]; Ds = 0;
X0 = Xs; X0_hat = Xs + [1;1;1;1];
T = 4; Ns = 200;
tspan = 0:T:Ns*T;
k_val = 1:1:length(tspan);
u = [2*sin(0.025*k_val)+0.15*cos(0.02*k_val);2*sin(0.02*k_val)-0.1*cos(0.025*k_val)];
U = u + Us;
n = 4; m=2;
X = zeros(n,Ns+1); 
X(:,1) = X0; 
load("Continuous_time_linear_perturbation_model_without_distubance.mat")
Phi = expm(A_mat*T); 
Gamma_u = (Phi-eye(size(Phi)))*pinv(A_mat)*B_mat;

%% Open Loop Simulation

for k=1:Ns
    % Plant simulation
    f_sys = @(t, X) System_Dynamics_210100059(X, params, U(:,k), 0); 
    [~,Y] = ode45(f_sys,[0,T],X(:,k)); 
    X(:,k+1) = (Y(end,:))';
end

[X_pred_OLE, SSE_OLE] = open_loop_estimator();    
[X_pred_PPO, SSE_PPO] = Pole_placement_Observer();
[X_pred_KF_st, X_est_KF_st, SSE_KF_st] = Kalman_Steady_Observer();
[X_pred_KF_var, X_est_KF_var, SSE_KF_var] = Kalman_Varying_Observer();   

X_est_PPO = X_pred_PPO; X_est_OLE = X_pred_OLE; % Estimation and prediction same in OLE and PPO

fprintf('\n\n  Metrics for the Observers:\n');
fprintf('----------------------------------------------------------------------------------------------------\n');
fprintf('%-20s %-25s %-18s %-15s %-10s\n', 'Metric', 'Open Loop Estimator', 'PP Observer', 'KF Steady', 'KF Varying');
fprintf('---------------------------------------------------------------------------------------------------\n');

for i = 1:4
    fprintf('SSE error in height %d     %-20.4f %-18.4f %-18.4f %-15.4f\n', i, SSE_OLE(i), SSE_PPO(i), SSE_KF_st(i), SSE_KF_var(i));
end


%% Plotting

observers = {'OLE', 'PPO', 'KF Steady', 'KF Varying'};
X_est_all = {X_est_OLE, X_est_PPO, X_est_KF_st, X_est_KF_var};
X_pred_all = {X_pred_OLE, X_pred_PPO, X_pred_KF_st, X_pred_KF_var};
colors = {'b', 'r', 'g', 'm'};

figure;
set(gcf, 'WindowState', 'maximized'); 
sgtitle("Estimated States", "FontSize", 25, "FontWeight", "bold");

for i = 1:4
    subplot(2,2,i)
    hold on; grid on;
    set(gca,"FontSize",15)
    xlabel("Time (in s)","FontSize",20); 
    ylabel("Water level (in cm)","FontSize",20);
    title(sprintf("$h_%d$", i), "FontSize", 25, "Interpreter", "latex")
    plot(tspan,X(i,:), "LineWidth", 2, 'Color', [0 1 1], 'DisplayName', "True State");
    
    for j = 1:4 
        plot(tspan, X_est_all{j}(i, :), "LineWidth", 2, 'Color', colors{j}, 'DisplayName', observers{j});
    end
    
    legend("FontSize", 12, "interpreter", "latex")
end

figure;
set(gcf, 'WindowState', 'maximized'); 
sgtitle("Manipulated Inputs", "FontSize", 25, "FontWeight", "bold");
subplot(2,1,1)
stairs(tspan,U(1,:),"LineWidth",2,"Color",[0 0 1])
xlabel("Time (in s)","FontSize",20); ylabel("Voltage level (in Volts)","FontSize",20); title("Input v1","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
subplot(2,1,2)
stairs(tspan,U(2,:),"LineWidth",2,"Color",[1 1 0])
xlabel("Time (in s)","FontSize",20); ylabel("Voltage level (in Volts)","FontSize",20); title("Input v1","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);

figure;
set(gcf, 'WindowState', 'maximized'); 
sgtitle("Predicted States", "FontSize", 25, "FontWeight", "bold");

for i = 1:4
    subplot(2,2,i)
    hold on; grid on;
    set(gca,"FontSize",15)
    xlabel("Time (in s)","FontSize",20); 
    ylabel("Water level (in cm)","FontSize",20);
    title(sprintf("$h_%d$", i), "FontSize", 25, "Interpreter", "latex")
    plot(tspan,X(i,:), "LineWidth", 2, 'Color', [0 1 1], 'DisplayName', "True State");
    
    for j = 1:4 
        plot(tspan, X_pred_all{j}(i, :), "LineWidth", 2, 'Color', colors{j}, 'DisplayName', observers{j});
    end
    
    legend("FontSize", 12, "interpreter", "latex")
end

