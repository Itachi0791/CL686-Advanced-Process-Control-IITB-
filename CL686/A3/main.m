clear;clc;close all;
T = 4; % Sampling Time (in s)
Ns = 200; % Number of Simulation steps
Us = [3;3]; Xs = [12.263;12.7831;1.6339;1.409];
U_L = [0;0]; U_H = [5;5];
R = Xs; X_L = [10;10;1;1]; X_H = [16;16;8;8]; 
[X_mpc,U_mpc,SSE_mpc,SSMV_mpc] = MPC();
[X_ppc,U_ppc,SSE_ppc,SSMV_ppc] = Pole_Placement_Control();
[X_lqr_var,U_lqr_var,SSE_lqr_var,SSMV_lqr_var] = LQOC_varying();
[X_lqr_st,U_lqr_st,SSE_lqr_st,SSMV_lqr_st] = LQOC_steady();

fprintf('\n\n  Metrics for the Controllers:\n');
fprintf('----------------------------------------------------------------------\n');
fprintf('%-28s %-10s %-8s %-10s %-10s\n', 'Metric', 'MPC', 'PPC', 'LQR Var', 'LQR St');
fprintf('----------------------------------------------------------------------\n');

for i = 1:length(SSE_mpc)
    fprintf('SSE error in height %d      %-10.4f %-10.4f %-10.4f %-10.4f\n', i, SSE_mpc(i), SSE_ppc(i), SSE_lqr_var(i), SSE_lqr_st(i));
end

for i = 1:length(SSMV_mpc)
    fprintf('SSMV error in voltage %d     %-10.4f %-10.4f %-10.4f %-10.4f\n', i, SSMV_mpc(i), SSMV_ppc(i), SSMV_lqr_var(i), SSMV_lqr_st(i));
end

%% Plotting

time = 0:T:(Ns-1)*T;
controllers = {'MPC', 'Pole Placement', 'LQR Time Varying', 'LQR Steady'};
X_all = {X_mpc, X_ppc, X_lqr_var, X_lqr_st};
U_all = {U_mpc, U_ppc, U_lqr_var, U_lqr_st};
colors = {'b', 'r', 'g', 'm'};

figure;
set(gcf, 'WindowState', 'maximized'); 
sgtitle("Controlled States", "FontSize", 25, "FontWeight", "bold");

for i = 1:4
    subplot(2,2,i)
    hold on; grid on;
    set(gca,"FontSize",15)
    xlabel("Time (in s)","FontSize",20); 
    ylabel("Water level (in cm)","FontSize",20);
    title(sprintf("$h_%d$", i), "FontSize", 25, "Interpreter", "latex")
    
    for j = 1:4 
        plot(time, X_all{j}(i, :), "LineWidth", 2, 'Color', colors{j}, 'DisplayName', controllers{j});
    end
    
    plot(time, R(i) * ones(size(time)), "y--", "LineWidth", 2, 'DisplayName', 'Setpoint'); 
    plot(time, X_L(i) * ones(size(time)), "k-.", "LineWidth", 2, 'DisplayName', 'Lower Bound');
    plot(time, X_H(i) * ones(size(time)), "k-.", "LineWidth", 2, 'DisplayName', 'Upper Bound');
    
    legend("FontSize", 12, "interpreter", "latex")
end

figure;
set(gcf, 'WindowState', 'maximized'); 
sgtitle("Manipulated Inputs", "FontSize", 25, "FontWeight", "bold");

for i = 1:2 
    subplot(2,1,i)
    hold on; grid on;
    set(gca,"FontSize",15)
    xlabel("Time (in s)","FontSize",20); 
    ylabel("Voltage level (in Volts)","FontSize",20);
    title(sprintf("$v_%d$", i), "FontSize", 25, "Interpreter", "latex")
    
    for j = 1:4 
        stairs(time(1:end-1), U_all{j}(i,:), "LineWidth", 2, 'Color', colors{j}, 'DisplayName', controllers{j});
    end
    
    plot(time(1:end-1), U_L(i) * ones(length(time)-1,1), "c-.", "LineWidth", 2, 'DisplayName', 'Lower Bound');
    plot(time(1:end-1), U_H(i) * ones(length(time)-1,1), "y-.", "LineWidth", 2, 'DisplayName', 'Upper Bound');
    
    legend("FontSize", 12, "interpreter", "latex")
end