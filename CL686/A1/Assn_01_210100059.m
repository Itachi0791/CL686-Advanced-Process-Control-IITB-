clear; clc; close all;
%% Physical Parameters for the Quadruple Tank Problem
params.A1 = 28; params.A2 = 32; params.A3 = 28; params.A4 = 32;
params.a1 = 0.071; params.a2 = 0.057; params.a3 = 0.071; params.a4 = 0.057;
params.g = 981;
params.k1 = 3.14; params.k2 = 3.29;
params.gamma1 = 0.43; params.gamma2 = 0.34; params.gamma3 = 0.4;

%% Simulation Variables
T = 4; % Sampling Time
Nstep = 300; Nprbs = 2000; % Number of Samples for step and prbs input respectively
% Steady state Conditions for the above parameters
Us = [3.15;3.15]; Ds = 2; Xs = [14.3;16.84;5.9;7.33];

load("Continuous_time_linear_model_with_disturbance.mat")
Phi = expm(A_mat*T); 
Gamma_u = (Phi-eye(size(Phi)))*pinv(A_mat)*B_mat;
Gamma_d = (Phi-eye(size(Phi)))*pinv(A_mat)*H_mat ;

% Initialize all matrices
Xstep = zeros(4,Nstep); Xstep_L = zeros(4,Nstep);
Xprbs = zeros(4,Nprbs); Xprbs_L = zeros(4,Nprbs);
Ustep = zeros(2,Nstep-1); Uprbs = zeros(2,Nprbs-1);
Dstep = zeros(Nstep-1,1); Dprbs = zeros(Nprbs-1,1);

% Set Initial Conditions
Xstep(:,1) = Xs; Xstep_L(:,1) = Xs;
Xprbs(:,1) = Xs; Xprbs_L(:,1) = Xs;

%% Step Input Simulation
x = zeros(4,Nstep); % To store perturbations for linear model
for k = 1 : Nstep-1
   if k < 50
       u_k = [0;0];
   else
       u_k = [0.5;0.5];
   end
   d_k = randn;
   Ustep(:,k) = Us + u_k; Dstep(k) = Ds + d_k;
   f = @(t, X) System_Dynamics_210100059(X, params, Ustep(:,k), Dstep(k)); 
   [~,Y] = ode45(f,[0,T],Xstep(:,k)); % One step non linear integration
   Xstep(:,k+1) = (Y(end,:))' ;
   x(:,k+1) = Phi*x(:,k) + Gamma_u*u_k + Gamma_d*d_k; % Linear Perturbation Model
   Xstep_L(:,k+1) = Xs + x(:,k+1); 
end
%% PRBS Input Simulation
u = (idinput([Nprbs,2],'prbs',[0,0.05],[-0.25,0.25]))'; % generating control pertubations in one go
x = zeros(4,Nprbs); % To store perturbations for linear model
for k = 1 : Nprbs-1
   d_k = randn;
   Uprbs(:,k) = Us + u(:,k); Dprbs(k) = Ds + d_k;
   f = @(t, X) System_Dynamics_210100059(X, params, Uprbs(:,k), Dprbs(k)); 
   [~,Y] = ode45(f,[0,T],Xprbs(:,k));  % One step non linear integration
   Xprbs(:,k+1) = (Y(end,:))' ;
   x(:,k+1) = Phi*x(:,k) + Gamma_u*u(:,k) + Gamma_d*d_k; % Linear Perturbation Model
   Xprbs_L(:,k+1) = Xs + x(:,k+1); 
end

%% Plotting all graphs
close all;
%Step input
time = 0:T:(Nstep-1)*T;
figure;
hold on
grid on
plot(time, Xstep(1, :), "LineWidth", 3, "Color", [0, 0, 1]);  
plot(time, Xstep(2, :), "LineWidth", 3, "Color", [1, 0, 0]);          
plot(time, Xstep(3, :), "LineWidth", 3, "Color", [0, 1, 0]);          
plot(time, Xstep(4, :), "LineWidth", 3, "Color", [0.7, 0.3, 0.8]);    
plot(time, Xstep_L(1, :), "LineWidth", 3, "LineStyle", ":", "Color", [1, 1, 0]);  
plot(time, Xstep_L(2, :), "LineWidth", 3, "LineStyle", ":", "Color", [0, 1, 1]); 
plot(time, Xstep_L(3, :), "LineWidth", 3, "LineStyle", ":", "Color", [0, 0, 0]);  
plot(time, Xstep_L(4, :), "LineWidth", 3, "LineStyle", ":", "Color", [1, 0.5, 0]);
legend({"h_1","h_2","h_3","h_4","h_1(Linear)","h_2(Linear)","h_3(Linear)","h_4(Linear)"},"NumColumns",2,"FontSize",17,"location","best")
xlabel("Time (in s)","FontSize",20); ylabel("Water level (in cm)","FontSize",20); title("Actual v/s Linear Model dynamics for Step Input","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
figure;
grid on
hold on
plot(time,Xstep(1,:)-Xstep_L(1,:),"LineWidth",2);
plot(time,Xstep(2,:)-Xstep_L(2,:),"LineWidth",2);
plot(time,Xstep(3,:)-Xstep_L(3,:),"LineWidth",2);
plot(time,Xstep(4,:)-Xstep_L(4,:),"LineWidth",2);
legend("$X_1-X_{L1}$","$X_2-X_{L2}$","$X_3-X_{L3}$","$X_4-X_{L4}$","FontSize",17,"NumColumns",2,"interpreter","latex")
xlabel("Time (in s)","FontSize",20); ylabel("Water level error (in cm)","FontSize",20); title("Error between Actual and Linear Model for Step Input","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
figure;
grid on
hold on
xlabel("Time (in s)","FontSize",20); ylabel("Voltage level (in Volts)","FontSize",20); title("Manipulated Inputs for Step input case","FontSize",20,"FontWeight","bold")
stairs(time(1:end-1),Ustep(1,:),"LineWidth",3,"Color",[0 0 1])
stairs(time(1:end-1),Ustep(2,:),"LineWidth",3,"LineStyle",":","Color",[1 1 0])
legend("v_1","v_2","FontSize",17)
set(gca, 'FontSize', 20);
figure;
hold on
grid on
stairs(time(1:end-1),Dstep(:),"LineWidth",2)
xlabel("Time (in s)","FontSize",20); ylabel("Disturbance flow (cm^3/s)","FontSize",20); title("Disturbance inputs for Step Input case","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);

% PRBS input
time = 0:T:(Nprbs-1)*T;
figure;
hold on
grid on
plot(time, Xprbs(1, :), "LineWidth", 3, "Color", [0, 0, 1]);  
plot(time, Xprbs(2, :), "LineWidth", 3, "Color", [1, 0, 0]);          
plot(time, Xprbs(3, :), "LineWidth", 3, "Color", [0, 1, 0]);          
plot(time, Xprbs(4, :), "LineWidth", 3, "Color", [0.7, 0.3, 0.8]);    
plot(time, Xprbs_L(1, :), "LineWidth", 3, "LineStyle", ":", "Color", [1, 1, 0]);  
plot(time, Xprbs_L(2, :), "LineWidth", 3, "LineStyle", ":", "Color", [0, 1, 1]); 
plot(time, Xprbs_L(3, :), "LineWidth", 3, "LineStyle", ":", "Color", [0, 0, 0]);  
plot(time, Xprbs_L(4, :), "LineWidth", 3, "LineStyle", ":", "Color", [1, 0.5, 0]);
legend({"h_1","h_2","h_3","h_4","h_1(Linear)","h_2(Linear)","h_3(Linear)","h_4(Linear)"},"NumColumns",2,"FontSize",17,"location","best")
xlabel("Time (in s)","FontSize",20); ylabel("Water level (in cm)","FontSize",20); title("Actual v/s Linear Model dynamics for PRBS Input","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
figure;
hold on
grid on
plot(time,Xprbs(1,:)-Xprbs_L(1,:),"LineWidth",2);
plot(time,Xprbs(2,:)-Xprbs_L(2,:),"LineWidth",2);
plot(time,Xprbs(3,:)-Xprbs_L(3,:),"LineWidth",2);
plot(time,Xprbs(4,:)-Xprbs_L(4,:),"LineWidth",2);
legend("$X_1-X_{L1}$","$X_2-X_{L2}$","$X_3-X_{L3}$","$X_4-X_{L4}$","FontSize",17,"NumColumns",2,"interpreter","latex")
xlabel("Time (in s)","FontSize",20); ylabel("Water level error (in cm)","FontSize",20); title("Error between Actual and Linear Model for PRBS Input","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
figure;
hold on
grid on
stairs(time(1:end-1),Uprbs(1,:),"LineWidth",2,"Color",[0 0 1])
stairs(time(1:end-1),Uprbs(2,:),"LineWidth",2,"Color",[1 1 0])
legend("v_1","v_2","FontSize",17)
xlabel("Time (in s)","FontSize",20); ylabel("Voltage level (in Volts)","FontSize",20); title("Manipulated Inputs for PRBS input case","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);
figure;
hold on
grid on
stairs(time(1:end-1),Dprbs(:),"LineWidth",2)
xlabel("Time (in s)","FontSize",20); ylabel("Disturbance flow (cm^3/s)","FontSize",20); title("Disturbance inputs for PRBS Input case","FontSize",20,"FontWeight","bold")
set(gca, 'FontSize', 20);