clc;clear;close all;
N_l = 3;
params.N_l = N_l;
n = 2*params.N_l + 4; m = params.N_l-1;
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
Phi0 = 0.01*ones(N_l-1,1);
Phi_dot0 = 0.01*ones(N_l-1,1);
px0 = 0;
py0 = 0.5;
theta0 = 0.5;
Z0 = [Phi0;theta0;px0;py0;Phi_dot0;0;0;0];
T = 0.1; tspan = 0:T:100;

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
% [t,Z] = ode45(@(t,Z) snake_dyn(Z,params,undulator(t,Z,params,Kp,Kd,Phi_dd_ref_fun,Phi_d_ref_fun,Phi_ref_fun)),tspan,Z0);
% Z_pred = zeros(size(Z)); Z_pred(1,:)=Z0; U = zeros(length(tspan)-1,m);
% for k = 1:length(tspan)-1
%     U(k,:) = undulator(t(k),Z(k,:)',params,Kp,Kd,Phi_dd_ref_fun,Phi_d_ref_fun,Phi_ref_fun);
%     Z_pred(k+1,:) = Z_pred(k,:) + T*snake_dyn(Z_pred(k,:)',params,U(k,:)')';
% end
% figure;
% plot(t(1:end-1),U)
% figure;
% plot(t,Z,"LineWidth",2);hold on; plot(t,Z_pred,"LineStyle","--","LineWidth",1.5)
%% MPC parameters and constraints

T = 0.01; % Sampling Time (in s)
p = 10; % Prediction Horizon   
q = 20; % Control Horizon
Us = zeros(N_l-1,1); px_s=0.1;py_s=0;theta_s=0; % State Setpoint
Xs = [zeros(N_l-1,1);theta_s;px_s;py_s;zeros(N_l+2,1)];  % Control Setpoint
n = size(Xs,1); m = size(Us,1); % No. of states and control inputs respectively 
wx = eye(n); % State Weights matrix
wu = 0.0000001*eye(m); % Control weights matrix
U_L = -1*ones(m,1); U_H = 1*ones(m,1); % Control input constraints
delU_L = -0.1*ones(m,1); delU_H = 0.1*ones(m,1); % Control change constraints
X_L = -100*ones(n,1); X_H = 100*ones(n,1); % State constraints
Phi_mat = @(A_mat,T) expm(A_mat*T); 

%% MPC Constant matrices

Wx = kron(eye(p),wx); % using kronecker product we can stack wx in block diagonal matrix p times
Wu = kron(eye(p),wu); % same for wu
Lambda = zeros(m*p,m*p); Lambda_0 = [eye(m);zeros((p-1)*m,m)]; % Matrices in change of control input inequalities
Lambda(1:m,1:m) = eye(m);
for k = 2:p
    Lambda((k-1)*m + 1:k*m,(k-2)*m + 1:(k-1)*m) = -eye(m);
    Lambda((k-1)*m + 1:k*m,(k-1)*m + 1:k*m) = eye(m);
end
Utilde_L = repmat(U_L,p,1); Utilde_H = repmat(U_H,p,1); % Control constraints stacked
Xtilde_L = repmat(X_L,p,1); Xtilde_H = repmat(X_H,p,1); % State constraints stacked
DeltaU_L = repmat(delU_L,p,1); DeltaU_H = repmat(delU_H,p,1); % Control change constraints stacked
Zs = repmat(Xs,p,1);
Aeq = zeros(m*(p-q),m*p); 
for i = 1:(p-q)
    Aeq((i-1)*m+1:i*m,(q-1)*m+1:q*m) = eye(m);
    Aeq((i-1)*m+1:i*m,(q+i-1)*m+1:(q+i)*m) = -eye(m);
end    
beq = zeros(m*(p-q),1);
A_inequality = @(Su)[eye(m*p);Lambda;Su;-eye(m*p);-Lambda;-Su];
b_inequality = @(x_k,u_kminus1,Del_k,Sx,Sdel) [Utilde_H;DeltaU_H+Lambda_0*u_kminus1;Xtilde_H-Sx*x_k-Sdel*Del_k;...
                                         -Utilde_L;-DeltaU_L-Lambda_0*u_kminus1;-Xtilde_L+Sx*x_k+Sdel*Del_k];

%% MPC simulation

Ns = 1000;
tspan = 0:T:Ns*T;
linearized(N_l) % Generates functions for A_matrix and B_matrix 
B_mat = B_matrix(params.m); % Constant
Z = zeros(n,Ns+1); Z(:,1)=Z0; U = zeros(m,length(tspan)-1);
% First step simulation with U(0) = 0
[~,Y] = ode45(@(t,Z) snake_dyn(Z,params,U(:,1)),[0,T],Z(:,1)); 
Z(:,2) = (Y(end,:))';
Ufk_prev = zeros(m*p,1);
options = optimoptions('quadprog','Algorithm','active-set','Display','off');
fprintf("Starting MPC Computation : \n")
for k = 2:Ns
    % To display progress
    if mod(k, round(Ns/5)) == 0
        fprintf('%d / %d iterations done\n', k, Ns);
    end
    A_mat= A_matrix(Z(:,k),params.m,params.c_p,params.c_n,params.c_t,params.lambda1,params.lambda2);
    Phi_mat = expm(A_mat*T);
    Gamma_mat = Gamma_matrix(A_mat, B_mat, T, n, m);
    del_k = T*snake_dyn(Z(:,k),params,U(:,k-1)) + Z(:,k) -Phi_mat*Z(:,k)-Gamma_mat*U(:,k-1);
    Del_k = repmat(del_k,p,1);
    [Sx,Su,Sdel,H,F] = MPC_matrices(Phi_mat,Gamma_mat,p,Wx,Wu,Z(:,k),Del_k,Zs);
    A_ineq = A_inequality(Su);
    b_ineq = b_inequality(Z(:,k),U(:,k-1),Del_k,Sx,Sdel);
    Ufk = quadprog(H,F,[],[],[],[],[],[], Ufk_prev,options);
    U(:,k) = Lambda_0'*Ufk;
    [~,Y] = ode45(@(t,Z) snake_dyn(Z,params,U(:,k)),[0,T],Z(:,k)); 
    Z(:,k+1) = (Y(end,:))';
end
plot(tspan(1:end-1),U)

%hold on
%plot(t,Z(:,N_l+1))


%% Animation
Z = Z';
Phi = Z(:,1:N_l-1);theta = Z(:,N_l); px = Z(:,N_l+1); py = Z(:,N_l+2);
t_cap = [cos(theta), sin(theta)]; n_cap = [-sin(theta), cos(theta)];
heads = cell(1,N_l); tails = cell(1,N_l); head = [px,py]; links = cell(1,N_l);
figure;hold on;grid on;
sgtitle("Snake Robot Undulation Control", "FontSize", 25, "FontWeight", "bold");
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

% %% Animation with realistic snake appearance
% Z = Z_pred;
% Phi = Z(:,1:N_l-1);
% theta = Z(:,N_l); 
% px = Z(:,N_l+1); 
% py = Z(:,N_l+2);
% t_cap = [cos(theta), sin(theta)]; 
% n_cap = [-sin(theta), cos(theta)];
% heads = cell(1,N_l); 
% tails = cell(1,N_l); 
% head = [px,py];
% 
% % Figure setup
% figure;
% hold on;
% set(gcf,'WindowState','maximized')
% xlim([-4, 4]); 
% ylim([-4, 4]);
% axis equal
% 
% % Calculate positions
% for i=1:N_l
%     heads{i} = head;
%     tails{i} = head - params.l*t_cap;
%     if i < N_l
%         head = tails{i} - Phi(:,i).*n_cap;
%     end
% end
% 
% % Initialize patches
% body_patch = patch('XData', [], 'YData', [], 'FaceColor', [0 0.7 0], 'EdgeColor', 'none');
% head_patch = patch('XData', [], 'YData', [], 'FaceColor', 'r', 'EdgeColor', 'none');
% 
% % Animation loop
% for j = 1:length(tspan)
%     % Get coordinates for this timestep
%     x_coords = zeros(1,N_l);
%     y_coords = zeros(1,N_l);
%     for i = 1:N_l
%         x_coords(i) = heads{i}(j,1);
%         y_coords(i) = heads{i}(j,2);
%     end
% 
%     % Create smooth curve
%     t = 1:N_l;
%     tt = 1:0.1:N_l;
%     xx = spline(t,x_coords,tt);
%     yy = spline(t,y_coords,tt);
% 
%     % Calculate width profile
%     base_width = params.l/3;
%     widths = linspace(base_width, base_width/3, length(tt));
% 
%     % Calculate normal vectors
%     dx = gradient(xx);
%     dy = gradient(yy);
%     norm = sqrt(dx.^2 + dy.^2);
%     normal_x = -dy./norm;
%     normal_y = dx./norm;
% 
%     % Create body polygon
%     upper_x = xx + widths.*normal_x;
%     upper_y = yy + widths.*normal_y;
%     lower_x = xx - widths.*normal_x;
%     lower_y = yy - widths.*normal_y;
% 
%     % Update body patch
%     set(body_patch, 'XData', [upper_x, fliplr(lower_x)], ...
%                     'YData', [upper_y, fliplr(lower_y)]);
% 
%     % Create head triangle
%     head_length = params.l;
%     head_width = base_width*5;
%     head_angle = theta(j);
% 
%     head_points_x = [
%         xx(1), ...                                           % Base center
%         xx(1) + head_length*cos(head_angle), ...            % Tip
%         xx(1) + head_width/2*normal_x(1), ...               % Right corner
%         xx(1) - head_width/2*normal_x(1)                    % Left corner
%     ];
% 
%     head_points_y = [
%         yy(1), ...                                          % Base center
%         yy(1) + head_length*sin(head_angle), ...           % Tip
%         yy(1) + head_width/2*normal_y(1), ...              % Right corner
%         yy(1) - head_width/2*normal_y(1)                   % Left corner
%     ];
% 
%     % Update head patch
%     set(head_patch, 'XData', head_points_x, ...
%                     'YData', head_points_y);
% 
%     drawnow;
%     pause(0.01*T);
% end
%% Functions

function Zdot = snake_dyn(Z,params,U)
    % Z = [Phi,theta,px,py,Phi_dot,theta_dot,px_dot,py_dot]^T
    N_l = params.N_l;
    m = params.m; 
    c_p = params.c_p;
    c_n = params.c_n; 
    c_t = params.c_t;
    lambda1 = params.lambda1; 
    lambda2 = params.lambda2;
    e_bar = params.e_bar;
    A = params.A; 
    D = params.D;

    Phi = Z(1:N_l-1); 
    theta = Z(N_l);
    px = Z(N_l + 1);
    py = Z(N_l + 2);
    Phi_dot = Z(N_l+3:2*N_l+1);  
    theta_dot = Z(2*N_l+2);
    px_dot = Z(2*N_l+3); 
    py_dot = Z(2*N_l+4);

    v_t = cos(theta)*px_dot + sin(theta)*py_dot;
    v_n = -sin(theta)*px_dot + cos(theta)*py_dot;

    Phi_ddot = (-c_n/m)*Phi_dot + (c_p/m)*v_t*A*D'*Phi+(1/m)*(D*D')*U;
    theta_ddot = (-lambda1)*theta_dot + (lambda2/(N_l-1))*v_t*e_bar'*Phi;
    v_t_dot = (-c_t/m)*v_t + (2*c_p/(N_l*m))*v_n*e_bar'*Phi - (c_p/(N_l*m))*Phi'*A*pinv(D)*Phi_dot;
    v_n_dot = (-c_n/m)*v_n + (2*c_p/(N_l*m))*v_t*e_bar'*Phi;
    px_ddot = v_t_dot*cos(theta) -v_n_dot*sin(theta) -v_t*sin(theta)*theta_dot -v_n*cos(theta)*theta_dot;
    py_ddot = v_t_dot*sin(theta) +v_n_dot*cos(theta) +v_t*cos(theta)*theta_dot -v_n*sin(theta)*theta_dot;

    Zdot = [Phi_dot;theta_dot;px_dot;py_dot;Phi_ddot;theta_ddot;px_ddot;py_ddot];
end

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

function Gamma = Gamma_matrix(A_matrix, B, T, n, m)
    % Convert B to a column vector to match ODE solver requirements
    B_vec = B(:); 
    
    % Define the ODE system
    ode_func = @(tau, x) reshape(expm(A_matrix * tau) * B, [], 1);
    
    % Solve the ODE
    [t, Y] = ode45(ode_func, [0, T], zeros(n * m, 1));
    
    % Extract the last value and reshape it back to matrix form
    Gamma = reshape(Y(end, :)', n, m);
end

function M = createSymbolicMatrix(name, m, n)
    % Create a symbolic matrix of size (m, n) with the given name
    M = sym(name, [m, n], 'real');
end
M1 = createSymbolicMatrix('A', 3, 4); % Creates a 3x4 symbolic matrix A
M2 = createSymbolicMatrix('B', 5, 2); % Creates a 5x2 symbolic matrix B
disp(M1);
disp(M2);

