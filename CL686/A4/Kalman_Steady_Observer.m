function [X_pred, X_est, SSE] = Kalman_Steady_Observer()    
    %clear; clc; close all;
        
    %% Physical Parameters for the Quadruple Tank Problem
    
    params.A1 = 28; params.A2 = 32; params.A3 = 28; params.A4 = 32;
    params.a1 = 0.071; params.a2 = 0.057; params.a3 = 0.071; params.a4 = 0.057;
    params.g = 981;
    params.k1 = 3.33; params.k2 = 3.35;
    params.gamma1 = 0.7; params.gamma2 = 0.6; params.gamma3 = 0.4;
    
    %% Simulation Parameters
    Xs = [12.263;12.7831;1.6339;1.409]; Us = [3;3]; Ds = 0;
    X0 = Xs; X0_est = Xs + [1;1;1;1];
    T = 4; Ns = 200;
    tspan = 0:T:Ns*T;
    k_val = 1:1:length(tspan);
    u = [2*sin(0.025*k_val)+0.15*cos(0.02*k_val);2*sin(0.02*k_val)-0.1*cos(0.025*k_val)];
    U = u + Us;
    n = 4; m = 2; n_op = 2;
    X = zeros(n,Ns+1); x_est = zeros(n,Ns+1);  
    X(:,1) = X0; x_est(:,1) = X0_est-Xs; x_pred = x_est;
    load("Continuous_time_linear_perturbation_model_without_distubance.mat")
    Phi = expm(A_mat*T); 
    Gamma = (Phi-eye(size(Phi)))*pinv(A_mat)*B_mat;

    %% Observer Design
    Q = 0.05^2*eye(n); R = 0.025^2*eye(n_op);    
    % Note, I am using discrete algebraic ricatti solver function, then
    % explicitly calculating L_inf in the form taught in class, since the
    % "kalman" function does not return L in the predictor-estimator form,
    % but only in the predictor form
    % Since estimate is a dual problem, we pass Phi transpose and C transpose as arguments,
    % it is solving the dual form of ricatti eqn in LQR
    [P_steady,~, ~] = idare(Phi', C_mat', Q, R, [], []);
    L = P_steady * C_mat' * pinv(C_mat * P_steady * C_mat' + R);
    
    %% Open Loop Simulation
    
    for k=1:Ns
        % Plant simulation
        f_sys = @(t, X) System_Dynamics_210100059(X, params, U(:,k), 0); 
        [~,Z] = ode45(f_sys,[0,T],X(:,k)); 
        X(:,k+1) = (Z(end,:))';
        y = C_mat*(X(:,k+1)-Xs); % Measurement
        % Prediction Step
        x_pred(:,k+1) = Phi*x_est(:,k) + Gamma*u(:,k);
        %Update Step
        x_est(:,k+1) = x_pred(:,k+1) + L*(y-C_mat*x_pred(:,k+1));
    end
    X_pred = x_pred + Xs;
    X_est = x_est + Xs;
    SSE = sum((X - X_est).^2, 2);
end