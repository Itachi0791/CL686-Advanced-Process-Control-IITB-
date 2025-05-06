clc; clear; close all;
%% Simulation Parameters

load("Continuous_time_linear_perturbation_model_without_distubance.mat")
Xs = [12.263;12.7831;1.6339;1.409]; Us = [3;3]; Ds = 0;
Ys = C_mat*Xs;
T = 4; Ns = 200;
tspan = 0:T:Ns*T;
y_sp = [ones(2,Ns/2) 10*ones(2,Ns/2+1)];
R = y_sp + Ys;
n = 4; m = 2; p = 2;
[X_inv,e_inv,U_inv,err_est_inv,Y_inv,SSEST_inv,SSESE_inv,SSEMV_inv] = innovation_bias();
[X_sa,e_sa,U_sa,err_est_sa,Y_sa,SSEST_sa,SSESE_sa,SSEMV_sa] = state_augmentation_input_bias();

fprintf('----------------------------------------------------------------------------------------------------\n');
fprintf('%-25s %-30s %-30s\n', 'Metric', 'LQOC with innovation Bias', 'LQOC with State Augmentation');
fprintf('---------------------------------------------------------------------------------------------------\n');

for i = 1:2
    fprintf('SSEST error in output %d         %-20.4f         %-20.4f\n', i, SSEST_inv(i), SSEST_sa(i));
end
for i = 1:4
    fprintf('SSESE error in height %d         %-20.4f         %-20.4f\n', i, SSESE_inv(i), SSESE_sa(i));
end
for i = 1:2
    fprintf('SSEMV error in input %d          %-20.4f         %-20.4f\n', i, SSEMV_inv(i), SSEMV_sa(i));
end


%% Plotting

figure;
set(gcf, 'WindowState', 'maximized'); 
sgtitle("States", "FontSize", 25, "FontWeight", "bold");

for i = 1:4
    subplot(2,2,i)
    hold on; grid on;
    set(gca,"FontSize",15)
    xlabel("Time (in s)","FontSize",20); 
    ylabel("Water level (in cm)","FontSize",20);
    title(sprintf("$h_%d$", i), "FontSize", 25, "Interpreter", "latex")
    plot(tspan,X_inv(i,:), "LineWidth", 2, 'Color', [1 0 0], 'DisplayName', "Innovation Bias");
    plot(tspan,X_sa(i,:), "LineWidth", 2, 'Color', [0 0 1], 'DisplayName', "State Augmentation");
    legend("FontSize", 12, "interpreter", "latex")
end

figure;
set(gcf, 'WindowState', 'maximized'); 
sgtitle("Estimation Error", "FontSize", 25, "FontWeight", "bold");

for i = 1:4
    subplot(2,2,i)
    hold on; grid on;
    set(gca,"FontSize",15)
    xlabel("Time (in s)","FontSize",20); 
    ylabel("Water level (in cm)","FontSize",20);
    title(sprintf("error$_%d$", i), "FontSize", 25, "Interpreter", "latex")
    plot(tspan,err_est_inv(i,:), "LineWidth", 2, 'Color', [1 0 0], 'DisplayName', "Innovation Bias");
    plot(tspan,err_est_sa(i,:), "LineWidth", 2, 'Color', [0 0 1], 'DisplayName', "State Augmentation");
    legend("FontSize", 12, "interpreter", "latex")
end

figure;
set(gcf, 'WindowState', 'maximized'); 
sgtitle("Innovations", "FontSize", 25, "FontWeight", "bold");
subplot(2,1,1)
hold on; grid on;
plot(tspan(1:end-1),e_inv(1,:),"LineWidth",2,"Color",[1 0 0], 'DisplayName', "Innovation Bias");
plot(tspan(1:end-1),e_sa(1,:),"LineWidth",2,"Color",[0 0 1], 'DisplayName', "State Augmentation");
xlabel("Time (in s)","FontSize",20); ylabel("Innovation level","FontSize",20); 
title("Innovation e1","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
legend("FontSize", 12, "interpreter", "latex")
subplot(2,1,2)
hold on; grid on;
plot(tspan(1:end-1),e_inv(2,:),"LineWidth",2,"Color",[1 0 0], 'DisplayName', "Innovation Bias");
plot(tspan(1:end-1),e_sa(2,:),"LineWidth",2,"Color",[0 0 1], 'DisplayName', "State Augmentation");
xlabel("Time (in s)","FontSize",20); ylabel("Innovation level","FontSize",20); 
title("Innovation e2","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
legend("FontSize", 12, "interpreter", "latex")

figure;
set(gcf, 'WindowState', 'maximized'); 
sgtitle("Manipulated Inputs", "FontSize", 25, "FontWeight", "bold");
subplot(2,1,1)
hold on; grid on;
stairs(tspan(1:end-1),U_inv(1,:),"LineWidth",2,"Color",[1 0 0], 'DisplayName', "Innovation Bias");
stairs(tspan(1:end-1),U_sa(1,:),"LineWidth",2,"Color",[0 0 1], 'DisplayName', "State Augmentation");
plot(tspan,Us(1)*ones(size(tspan)),"LineWidth",2,"Color",[0 1 0], 'DisplayName', "Steady state")
xlabel("Time (in s)","FontSize",20); ylabel("Voltage level (in Volts)","FontSize",20); 
title("Input v1","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
legend("FontSize", 12, "interpreter", "latex")
subplot(2,1,2)
hold on; grid on;
stairs(tspan(1:end-1),U_inv(2,:),"LineWidth",2,"Color",[1 0 0], 'DisplayName', "Innovation Bias");
stairs(tspan(1:end-1),U_sa(2,:),"LineWidth",2,"Color",[0 0 1], 'DisplayName', "State Augmentation");
plot(tspan,Us(2)*ones(size(tspan)),"LineWidth",2,"Color",[0 1 0], 'DisplayName', "Steady state")
xlabel("Time (in s)","FontSize",20); ylabel("Voltage level (in Volts)","FontSize",20); 
title("Input v2","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
legend("FontSize", 12, "interpreter", "latex")

figure;
set(gcf, 'WindowState', 'maximized'); 
sgtitle("Measurement Outputs", "FontSize", 25, "FontWeight", "bold");
subplot(2,1,1)
hold on; grid on;
plot(tspan,Y_inv(1,:),"LineWidth",2,"Color",[1 0 0], 'DisplayName', "Innovation Bias");
plot(tspan,Y_sa(1,:),"LineWidth",2,"Color",[0 0 1], 'DisplayName', "State Augmentation");
plot(tspan,R(1,:),"LineWidth",2,"Color",[0 1 0], 'DisplayName', "Setpoint");
xlabel("Time (in s)","FontSize",20); ylabel("Output level","FontSize",20); 
title("Output y1","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
legend("FontSize", 12, "interpreter", "latex")
subplot(2,1,2)
hold on; grid on;
plot(tspan,Y_inv(2,:),"LineWidth",2,"Color",[1 0 0], 'DisplayName', "Innovation Bias");
plot(tspan,Y_sa(2,:),"LineWidth",2,"Color",[0 0 1], 'DisplayName', "State Augmentation");
plot(tspan,R(2,:),"LineWidth",2,"Color",[0 1 0], 'DisplayName', "Setpoint");
xlabel("Time (in s)","FontSize",20); ylabel("Output level","FontSize",20); 
title("Output y2","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
legend("FontSize", 12, "interpreter", "latex")



