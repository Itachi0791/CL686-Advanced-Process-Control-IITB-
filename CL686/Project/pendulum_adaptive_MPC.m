clc; clear; close all;
g = 9.81; l = 1; alpha = g/l; b = 0;
params.alpha = alpha; params.b = b;
n = 2 ; m = 1;
Ns = 200;  % No. of time steps
T = 0.05; % Sampling Time
p = 100; % Prediction Horizon   
q = 20; % Control Horizon
tspan = 0:T:Ns*T;
x0 = [0;0];%[2*pi*rand-pi;6*rand-3]; 
wx = eye(n); % State Weights matrix
wu = 0.000*eye(m); % Control weights matrix
U_L = -10*ones(m,1); U_H = 10*ones(m,1); % Control input constraints
delU_L = -1.5*ones(m,1); delU_H = 1.5*ones(m,1); % Control change constraints
X_L = [-pi;-4] ; X_H = [pi;4]; % State constraints
Xs = [3*pi/4;0];%[2*pi*rand-pi;0]; % Setpoint

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

X = zeros(n,Ns+1); X(:,1)=x0; U = zeros(m,length(tspan)-1);
% First step simulation with U(0) = 0
[~,Y] = ode45(@(t,x) pendulum_dyn(x,params,U(:,1)),[0,T],X(:,1)); 
X(:,2) = (Y(end,:))';
Ufk_prev = zeros(m*p,1);
options = optimoptions('quadprog','Algorithm','active-set','Display','off');
fprintf("Starting MPC Computation : \n")
for k = 2:Ns
     % To display progress
    if mod(k, round(Ns/5)) == 0
        fprintf('%d / %d iterations done\n', k, Ns);
    end
    B_mat = B_matrix(X(:,k),params);
    A_mat= A_matrix(X(:,k),params);
    Phi_mat = Phi_matrix(A_mat,T);
    Gamma_mat = Gamma_matrix(A_mat, B_mat, T, n, m);
    del_k = T*pendulum_dyn(X(:,k),params,U(:,k-1)) + X(:,k) -Phi_mat*X(:,k)-Gamma_mat*U(:,k-1);
    Del_k = repmat(del_k,p,1);
    [Sx,Su,Sdel,H,F] = MPC_matrices(Phi_mat,Gamma_mat,p,Wx,Wu,X(:,k),Del_k,Zs);
    A_ineq = A_inequality(Su);
    b_ineq = b_inequality(X(:,k),U(:,k-1),Del_k,Sx,Sdel);
    Ufk = quadprog(H,F,A_ineq,b_ineq,Aeq,beq,[],[], Ufk_prev,options);
    U(:,k) = Lambda_0'*Ufk;
    [~,Y] = ode45(@(t,x) pendulum_dyn(x,params,U(:,k)),[0,T],X(:,k)); 
    X(:,k+1) = (Y(end,:))';
    Ufk_prev = Ufk;
end
fprintf("MPC Computation completed.\n")
%% Plots and Animation
close all;
figure;
sgtitle("Adaptive MPC - Controlled States", "FontSize", 25, "FontWeight", "bold");
subplot(2,1,1)
plot(tspan, X(1,:),LineWidth=2);
hold on;
plot(tspan,Xs(1)*ones(size(tspan)),"LineStyle","-.",LineWidth=2);
plot(tspan,X_H(1)*ones(size(tspan)),"k","LineStyle","--",LineWidth=2);
plot(tspan,X_L(1)*ones(size(tspan)),"k","LineStyle","--",LineWidth=2);
ylabel('$\theta$','Interpreter','latex','FontSize',25,'FontWeight','bold','Rotation',0)
legend('Actual', 'Reference', 'Constraints','FontSize',12)
grid on
set(gca,'fontsize',18)
subplot(2,1,2)
plot(tspan, X(2,:),LineWidth=2);
hold on;
plot(tspan, Xs(2)*ones(size(tspan)),"LineStyle","-.",LineWidth=2);
plot(tspan,X_H(2)*ones(size(tspan)),"k","LineStyle","--",LineWidth=2);
plot(tspan,X_L(2)*ones(size(tspan)),"k","LineStyle","--",LineWidth=2);
xlabel("Time [in s]",'FontSize',17)
ylabel('$\dot{\theta}$','Interpreter','latex','FontSize',25,'FontWeight','bold','Rotation',0)
legend('Actual', 'Reference', 'Constraints','FontSize',12)
grid on
set(gca,'fontsize',18)

figure;
stairs(tspan(1:end-1),U,LineWidth=2)
hold on
plot(tspan,U_H(1)*ones(size(tspan)),"k","LineStyle","--",LineWidth=2);
plot(tspan,U_L(1)*ones(size(tspan)),"k","LineStyle","--",LineWidth=2);
ylabel('U','FontSize',25,'FontWeight','bold','Rotation',0)
xlabel("Time [in s]",'FontSize',17)
legend('Actual','Constraints','FontSize',12)
title("Adaptive MPC - Manipulated Inputs", "FontSize", 25, "FontWeight", "bold");
grid on
set(gca,'fontsize',18)

% % Visualize and animate pendulum motion
% figure('Position', [100, 100, 700, 525]);
% hold on;
% 
% % Initialize pendulum graphics
% bob = plot(0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % Pendulum bob
% rod = line([0, 0], [0, 0], 'Color', 'k', 'LineWidth', 2); % Pendulum rod
% time_text = text(-l * 1.1, l * 1.1, 'Time: 0.00 s', 'FontSize', 18, 'Color', 'r'); % Time text
% 
% % Set axis limits and aspect ratio
% axis equal;
% set(gca,'FontSize',20)
% grid on;
% axis([-l, l, -l, l] * 1.2);
% xlabel('X axis (m)');
% ylabel('Y axis (m)');
% title('Adaptive MPC for Pendulum');
% theta = X(1,:);
% 
% for k = 1:length(theta)
%     % Calculate pendulum position
%     xpos = l * sin(theta(k));
%     ypos = -l * cos(theta(k));
% 
%     % Update graphics
%     set(bob, 'XData', xpos, 'YData', ypos); % Update bob position
%     set(rod, 'XData', [0, xpos], 'YData', [0, ypos]); % Update rod position
%     set(time_text, 'String', sprintf('Time : %.2f s', (k-1)*T)); % Update time text
% 
%     % % Capture the current frame
%     % frame = getframe(gcf); % Capture the figure as a frame
%     % writeVideo(video, frame); % Write the frame to the video
% 
%     % Pause to control animation speed
%     pause(0.1);
% end
% 
% hold off;


%% Functions

function A = A_matrix(x,params)
    A = [0 1; -params.alpha*cos(x(1)) -params.b];
end

function B = B_matrix(x,params)
    B = [0;1];
end

function Phi = Phi_matrix(A_matrix,T)
    Phi = expm(A_matrix*T);
end

function Gamma = Gamma_matrix(A_matrix, B, T, n, m)
    if det(A_matrix)>1e-10
        Gamma = (expm(A_matrix*T)-eye(n))*pinv(A_matrix)*B;
    else
        ode_func = @(tau, x) reshape(expm(A_matrix * tau) * B, [], 1);
        [~, Y] = ode45(ode_func, [0, T], zeros(n * m, 1));
        Gamma = reshape(Y(end, :)', n, m);
    end
end

function xdot = pendulum_dyn(x,params,u)
    xdot = [x(2);-params.b*x(2)-params.alpha*sin(x(1))+u];
end
