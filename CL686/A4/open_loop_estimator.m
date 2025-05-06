function [X_hat, SSE] = open_loop_estimator()    
    %clear; clc; close all;
        
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
    X = zeros(n,Ns+1); x_hat = zeros(n,Ns+1);
    X(:,1) = X0; x_hat(:,1) = X0_hat-Xs;
    load("Continuous_time_linear_perturbation_model_without_distubance.mat")
    Phi = expm(A_mat*T); 
    Gamma_u = (Phi-eye(size(Phi)))*pinv(A_mat)*B_mat;
    
    %% Open Loop Simulation
    
    for k=1:Ns
        % Plant simulation
        f_sys = @(t, X) System_Dynamics_210100059(X, params, U(:,k), 0); 
        [~,Y] = ode45(f_sys,[0,T],X(:,k)); 
        X(:,k+1) = (Y(end,:))';
        % Estimate Propagation
        x_hat(:,k+1) = Phi*x_hat(:,k) + Gamma_u*u(:,k);
    end
    
    X_hat = x_hat + Xs;
    SSE = sum((X - X_hat).^2, 2);
end