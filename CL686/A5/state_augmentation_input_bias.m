function [X,e,U,err_est,Y,SSEST,SSESE,SSEMV] = state_augmentation_input_bias()
    %clear; clc; close all;
    
    %% Physical Parameters for the Quadruple Tank Problem
    
    params.A1 = 28; params.A2 = 32; params.A3 = 28; params.A4 = 32;
    params.a1 = 0.071; params.a2 = 0.057; params.a3 = 0.071; params.a4 = 0.057;
    params.g = 981;
    params.k1 = 3.33; params.k2 = 3.35;
    params.gamma1 = 0.7; params.gamma2 = 0.6; params.gamma3 = 0.4;
    
    %% Simulation Parameters
    
    load("Continuous_time_linear_perturbation_model_without_distubance.mat")
    Xs = [12.263;12.7831;1.6339;1.409]; Us = [3;3]; Ds = 0;
    X0 = Xs; x0_est = [1;1;1;1];
    Ys = C_mat*Xs;
    T = 4; Ns = 200;
    tspan = 0:T:Ns*T;
    y_sp = [ones(2,Ns/2) 10*ones(2,Ns/2+1)];
    R = y_sp + Ys;
    n = 4; m = 2; p = 2;
    X = zeros(n,Ns+1); x_est = zeros(n,Ns+1); U = zeros(m,Ns); u = zeros(m,Ns);
    e = zeros(p,Ns); beta_est = zeros(p,Ns+1);
    X(:,1) = X0; x_est(:,1) = x0_est;
    Phi = expm(A_mat*T); 
    Gamma_u = (Phi-eye(size(Phi)))*pinv(A_mat)*B_mat;
    
    %% Controller Design
    
    Wx = 10*eye(n); Wu = eye(m);
    [G_inf,~,~] = dlqr(Phi,Gamma_u,Wx,Wu);
    
    %% Observer Design
    
    Q = 0.2*eye(n+p); R_mat = 5*eye(m);
    Phi_sa = [Phi Gamma_u; zeros(p,n) eye(p,p)];
    Gamma_sa = [Gamma_u;zeros(p,p)];
    C_sa = [C_mat zeros(p,p)];
    [L_T,~,~] = dlqr(Phi_sa', C_sa', Q, R_mat);
    L_inf_sa = L_T';

    %% State Augmentation (Input Bias)
    
    K_u = C_mat*pinv(eye(n)-Phi)*Gamma_u;
    K_u_inv = pinv(K_u); inv_IminPhi = pinv(eye(n)-Phi);
    x_sa = [x_est(:,1);beta_est(:,1)];
    
    for k=1:Ns
        % Innovation
        e(:,k) = C_mat*(X(:,k)-Xs-x_est(:,k));
        % Control
        u_sp = K_u_inv*y_sp(:,k)-beta_est(:,k);
        x_sp = inv_IminPhi*Gamma_u*K_u_inv*y_sp(:,k);
        u(:,k) = u_sp - G_inf*(x_est(:,k)-x_sp);
        % Estimation
        x_sa = Phi_sa*x_sa + Gamma_sa*u(:,k) + L_inf_sa*e(:,k);
        x_est(:,k+1) = x_sa(1:n); beta_est(:,k+1) = x_sa(n+1:n+p);
        % Plant simulation
        U(:,k) = u(:,k) + Us;
        f_sys = @(t, X) System_Dynamics_210100059(X, params, U(:,k), 0); 
        [~,Z] = ode45(f_sys,[0,T],X(:,k)); 
        X(:,k+1) = (Z(end,:))';
    end
    Y = C_mat*X;
    X_est = x_est + Xs;
    err_est = X - X_est; 
    SSEST = sum((Y-R).^2,2);
    SSESE = sum((X - X_est).^2, 2);
    SSEMV = sum((U-Us).^2,2);
end

