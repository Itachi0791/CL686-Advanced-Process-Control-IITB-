clc;clear;close all;
N_l = 10;
params.N_l = N_l;
n = 2*params.N_l + 4; % Number of state variables
m = params.N_l-1; % Number of control inputs
A = zeros(N_l-1,N_l); D = A;
for k = 1:N_l-1
    A(k,k:k+1) = [1 1];
    D(k,k:k+1) = [1 -1];
end
params.A = A; params.D = D;
params.c_t = 1; params.c_n = 3;
params.lambda1 = 0.5; params.lambda2 = 20;
params.m = 1; params.l = 0.3; 
params.e = ones(N_l,1); params.e_bar = ones(N_l-1,1);
params.c_p = (params.c_t-params.c_n)/(2*params.l);
Phi0 = zeros(N_l-1,1);
Phi_dot0 = zeros(N_l-1,1);
px0 = 0;
py0 = 0.5;
theta0 = 0;
Z0 = [Phi0;theta0;px0;py0;Phi_dot0;0;0;0];
T = 0.01; tspan = 0:T:50;

Kp = 20;
Kd = 5;
w = (-2*pi/3);
delta = (2*pi/9);
alpha = 0.05;
Phi_mean = 0;

Phi_ref = zeros(length(tspan),N_l-1);
Phi_d_ref = zeros(length(tspan),N_l-1);
Phi_dd_ref = zeros(length(tspan),N_l-1);

for j=1:length(tspan)
    t = tspan(j);
    for k=1:N_l-1
        Phi_ref(j,k) = Phi_mean + alpha * sin(w*t + (k - 1)*delta);
        Phi_d_ref(j,k) = alpha*w*cos(w*t + (k - 1)*delta);
        Phi_dd_ref(j,k) = -alpha*w^2*sin(w*t + (k - 1)*delta);
    end
end

Phi_dd_ref_fun = @(t)interp1(tspan.',Phi_dd_ref,t).';
Phi_d_ref_fun = @(t)interp1(tspan.',Phi_d_ref,t).';
Phi_ref_fun = @(t)interp1(tspan.',Phi_ref,t).';

%% Euler Comparison
[t,Z] = ode45(@(t,Z) snake_dyn(Z,params,undulator(t,Z,params,Kp,Kd,Phi_dd_ref_fun,Phi_d_ref_fun,Phi_ref_fun)),tspan,Z0);
Z_pred = zeros(size(Z)); Z_pred(1,:)=Z0; U = zeros(length(tspan)-1,m);
for k = 1:length(tspan)-1
    U(k,:) = undulator(t(k),Z(k,:)',params,Kp,Kd,Phi_dd_ref_fun,Phi_d_ref_fun,Phi_ref_fun);
    Z_pred(k+1,:) = Z_pred(k,:) + T*snake_dyn(Z_pred(k,:)',params,U(k,:)')';
end
figure;
plot(t(1:end-1),U)
figure;
plot(t,Z,"LineWidth",2);hold on; plot(t,Z_pred,"LineStyle","--","LineWidth",1.5)
%% MPC parameters and constraints

Ts = 0.08; % Sampling Time (in s)
Ns = 400; % Number of Time samples


Np = 20; % Prediction Horizon 
tspan_pred = 0:Ts:Ts*(Np-1);
tspan = 0:Ts:Ts*Ns;
U = zeros(length(tspan)-1,m);
Z = zeros(length(tspan),n);
Ufk_start= 0.01*ones(Np,m); 

for i=1:Ns
    disp(i)
    objective = @(Ufk)MPC_objective(Ts,Np,Z(i,:),params,Ufk);
    constraints = @(Ufk)MPC_constraints(Ts,Np,Z(i,:),params,Ufk);
    options = optimoptions('fmincon','Display','off','Algorithm','sqp','PlotFcn','optimplotfval','MaxFunctionEvaluations',1e4);
    Ufk = fmincon(objective,Ufk_start,[],[],[],[],-0.2*ones(Np,m),0.2*ones(Np,m),[],options);
    Ufk_start = Ufk;
    U(i,:) = Ufk(1,:);
    [~,Y] = ode45(@(t,Z) snake_dyn(Z,params,U(i,:)'),[0,Ts],Z(i,:)'); 
    Z(i+1,:) = (Y(end,:))';
end

figure;
plot(tspan_pred, Ufk);

%% Simulation with MPC


[t,Z_pred] = ode45(@(t,Z) snake_dyn_2(t,Z,params,tspan_pred,Ufk),tspan_pred,Z0);

figure;
plot(tspan_pred, Z_pred(:,2*N_l+3));
hold on
plot(tspan_pred, Z_pred(:,2*N_l+4));

%% Animation
Z_pred = Z;
Phi = Z_pred(:,1:N_l-1);theta = Z_pred(:,N_l); px = Z_pred(:,N_l+1); py = Z_pred(:,N_l+2);
t_cap = [cos(theta), sin(theta)]; n_cap = [-sin(theta), cos(theta)];

% Phi = Z(:,1:N_l-1);theta = Z(:,N_l); px = Z(:,N_l+1); py = Z(:,N_l+2);
% t_cap = [cos(theta), sin(theta)]; n_cap = [-sin(theta), cos(theta)];

heads = cell(1,N_l); tails = cell(1,N_l); head = [px,py]; links = cell(1,N_l);
figure;hold on;
set(gcf,'WindowState','maximized')
xlim([-4, 4]); ylim([-4, 4]);

for k=1:N_l
    heads{k} = head;
    tails{k} = head - params.l*t_cap;
    if k < N_l
        head = tails{k} - Phi(:,k).*n_cap;
        links{k} = plot([heads{k}(1,1) head(1,1)],[heads{k}(1,2) head(1,2)], 'b-', 'LineWidth', 4);
    end
    %links{i} = plot([heads{i}(1,1) tails{i}(1,1)],[heads{i}(1,2) tails{i}(1,2)], 'b-', 'LineWidth', 2);
end

for j = 1:length(tspan)
    for k = 1:N_l-1
        set(links{k},'Xdata',[heads{k}(j,1) heads{k+1}(j,1)],'Ydata',[heads{k}(j,2) heads{k+1}(j,2)]);
        %set(links{i},'Xdata',[heads{i}(j,1) tails{i}(j,1)],'Ydata',[heads{i}(j,2) tails{i}(j,2)]);
    end
    pause(0.01);
end

%% Functions

function U = undulator(t,Z,params,Kp,Kd,Phi_dd_ref_fun, Phi_d_ref_fun, Phi_ref_fun)
    A = params.A;
    D = params.D;
    cn = params.c_n;
    cp = params.c_p;
    N_l = params.N_l;
    m = params.m;
    theta = Z(N_l);
    Phi = Z(1:N_l-1);
    Phi_dot = Z(N_l+3:2*N_l+1);
    px_dot = Z(2*N_l+3); py_dot = Z(2*N_l+4);
    vt = cos(theta)*px_dot + sin(theta)*py_dot;

    U_bar = Phi_dd_ref_fun(t) + Kp * (Phi_ref_fun(t) - Z(1:params.N_l-1)) + Kd * (Phi_d_ref_fun(t) - Z(params.N_l+3: 2*params.N_l+1));
    U = m * inv(D*D') * (U_bar + (cn/m)*Phi_dot - (cp/m)*vt*A*D'*Phi);
end



