%% Deriving the equations of motion of a One-link manipulator with a mobile base using Symbolic  Math Toolbox


% In this code, the equations governing the dynamics of the
% loco-manipulator system (one link manipulator with a mobile base) are
% derived using the Symbolic Math Toolbox

% Link 1 = Mobile Base
% Link 2 = Manipulator Link

%% Create symbols

syms m1 l1 lc1 m2 l2 lc2 g  % Mass and Length parameters
syms posx(t) posy(t) posx_d(t) posy_d(t) theta1(t) theta1_d(t) theta2(t) theta2_d(t) % Placeholder state variables
syms t % Time
syms x0 x0_d x0_dd y0 y0_d y0_dd q1 q1_d q1_dd q2 q2_d q2_dd % State variables and their derivatives as a function of time
syms Fx Fy tau1 tau2 % Force and Torque actuator inputs (in case mobile has two independent thrusts along both x and y directions)
%syms Ft tau1 tau2

%% Create variables

I1 = (m1 * l1^2)/12; % Inertia of link 1 about COM
I2 = (m2 * l2^2)/12; % Inertia of link 2 about COM
w1 = theta1_d(t); % Angular velocity of link 1
w2 = theta1_d(t) + theta2_d(t); % Angular velocity of link 2

pos2x_d = posx_d(t) - l1 * w1 * sin(theta1(t)); % x-component of velocity of the joint between link 1 and link 2
pos2y_d = posy_d(t) + l1 * w1 * cos(theta1(t)); % y-component of velocity of the joint between link 1 and link 2

K1 = 0.5 * I1 * w1.^2 +  0.5 * m1 * ((lc1 * w1 * cos(theta1(t)) + posy_d(t))^2 + (-lc1 * w1 * sin(theta1(t)) + posx_d(t))^2); % KE of link 1
K2 = 0.5 * I2 * w2.^2 +  0.5 * m2 * ((lc2 * w2 * cos(theta1(t) + theta2(t)) + pos2y_d)^2 + (-lc2 * w2 * sin(theta1(t) + theta2(t)) + pos2x_d)^2); % KE of link 2
K = K1 + K2; % Total Kinetic energy

P1 = m1 * g * (posy(t) + lc1 * sin(theta1(t))); % PE of link 1
P2 = m2 * g * (posy(t) + (l1 * sin(theta1(t)) + lc2 * sin(theta1(t)+theta2(t)))); % PE of link 2
P = P1 + P2; % Total potential energy

L = K - P; % Lagrangian

%% Euler-Lagrange equation

%Fx = Ft * cos(theta1(t)); % x-component of thrust on the base
%Fy = Ft * sin(theta1(t)); % y-component of thrust on the base

Eqn1 = diff(diff(L,theta1_d(t)),t) - diff(L,theta1(t)) - tau1 == 0; % For link 1
Eqn2 = diff(diff(L,theta2_d(t)),t) - diff(L,theta2(t)) - tau2 == 0; % For link 2
Eqn3 = diff(diff(L,posx_d(t)),t) - diff(L,posx(t)) - Fx == 0; % For x-component of mobile base
Eqn4 = diff(diff(L,posy_d(t)),t) - diff(L,posy(t)) - Fy == 0; % For y-component of mobile base

% Substitute q1_dd, q2_dd, x0_dd, y0_dd into the equations
Eqn1 = subs(Eqn1, [diff(theta1_d(t),t),diff(theta2_d(t),t),diff(posx_d(t),t),diff(posy_d(t),t)], [q1_dd,q2_dd,x0_dd,y0_dd]);
Eqn2 = subs(Eqn2, [diff(theta1_d(t),t),diff(theta2_d(t),t),diff(posx_d(t),t),diff(posy_d(t),t)], [q1_dd,q2_dd,x0_dd,y0_dd]);
Eqn3 = subs(Eqn3, [diff(theta1_d(t),t),diff(theta2_d(t),t),diff(posx_d(t),t),diff(posy_d(t),t)], [q1_dd,q2_dd,x0_dd,y0_dd]);
Eqn4 = subs(Eqn4, [diff(theta1_d(t),t),diff(theta2_d(t),t),diff(posx_d(t),t),diff(posy_d(t),t)], [q1_dd,q2_dd,x0_dd,y0_dd]);

% Substitute q1_d, q2_d, x0_d, y0_d into the equations
Eqn1 = subs(Eqn1, [diff(theta1(t),t),diff(theta2(t),t),diff(posx(t),t),diff(posy(t),t)], [q1_d,q2_d,x0_d,y0_d]);
Eqn2 = subs(Eqn2, [diff(theta1(t),t),diff(theta2(t),t),diff(posx(t),t),diff(posy(t),t)], [q1_d,q2_d,x0_d,y0_d]);
Eqn3 = subs(Eqn3, [diff(theta1(t),t),diff(theta2(t),t),diff(posx(t),t),diff(posy(t),t)], [q1_d,q2_d,x0_d,y0_d]);
Eqn4 = subs(Eqn4, [diff(theta1(t),t),diff(theta2(t),t),diff(posx(t),t),diff(posy(t),t)], [q1_d,q2_d,x0_d,y0_d]);

Eqn1 = subs(Eqn1, [theta1_d(t),theta2_d(t),posx_d(t),posy_d(t)], [q1_d,q2_d,x0_d,y0_d]);
Eqn2 = subs(Eqn2, [theta1_d(t),theta2_d(t),posx_d(t),posy_d(t)], [q1_d,q2_d,x0_d,y0_d]);
Eqn3 = subs(Eqn3, [theta1_d(t),theta2_d(t),posx_d(t),posy_d(t)], [q1_d,q2_d,x0_d,y0_d]);
Eqn4 = subs(Eqn4, [theta1_d(t),theta2_d(t),posx_d(t),posy_d(t)], [q1_d,q2_d,x0_d,y0_d]);

% Substitute q1,q2,x0,y0 into the equations
Eqn1 = subs(Eqn1,[theta1(t),theta2(t),posx(t),posy(t)],[q1,q2,x0,y0]); 
Eqn2 = subs(Eqn2,[theta1(t),theta2(t),posx(t),posy(t)],[q1,q2,x0,y0]);
Eqn3 = subs(Eqn3,[theta1(t),theta2(t),posx(t),posy(t)],[q1,q2,x0,y0]);
Eqn4 = subs(Eqn4,[theta1(t),theta2(t),posx(t),posy(t)],[q1,q2,x0,y0]);

% Solve for q1_ddot, q2_ddot, x0_dd, y0_dd
[q1_dd, q2_dd, x0_dd, y0_dd] = solve([Eqn1, Eqn2, Eqn3, Eqn4], [q1_dd, q2_dd, x0_dd, y0_dd]);

q1_dd = simplify(q1_dd);
q2_dd = simplify(q2_dd);
x0_dd = simplify(x0_dd);
y0_dd = simplify(y0_dd);

% Expression for Power due to the actuator inputs

Pow = tau1 * q1_d + tau2 * q2_d + Fx * x0_d + Fy * y0_d;

% Substituting the placeholder state variables in the expressions for Kinetic Energy, Potential Energy and Power
K = subs(K,[theta1(t), theta1_d(t), theta2(t), theta2_d(t), posx(t), posx_d(t), posy(t), posy_d(t)],[q1, q1_d, q2, q2_d, x0, x0_d, y0, y0_d]);
P = subs(P,[theta1(t), theta1_d(t), theta2(t), theta2_d(t), posx(t), posx_d(t), posy(t), posy_d(t)],[q1, q1_d, q2, q2_d, x0, x0_d, y0, y0_d]);
Pow = subs(Pow,[theta1(t), theta1_d(t), theta2(t), theta2_d(t), posx(t), posx_d(t), posy(t), posy_d(t)],[q1, q1_d, q2, q2_d, x0, x0_d, y0, y0_d]);

% State derivative vector
dx = [q1_d; q2_d; x0_d; y0_d; q1_dd; q2_dd; x0_dd; y0_dd; Pow];

%% Create Matrices for Linearized System

A_mat = jacobian([q1_d; q2_d; x0_d; y0_d; q1_dd; q2_dd; x0_dd; y0_dd],[q1;q2;x0;y0;q1_d;q2_d;x0_d;y0_d]);
B_mat = jacobian([q1_d; q2_d; x0_d; y0_d; q1_dd; q2_dd; x0_dd; y0_dd],[tau1;tau2;Fx;Fy]);

%% Create function for state evolution

x_vec = [q1;q2;x0;y0;q1_d;q2_d;x0_d;y0_d];
u_vec = [tau1;tau2;Fx;Fy];

% Create function dx.m to use in ode
matlabFunction(dx,"File","dx","Vars",{t,g,l1,lc1,l2,lc2,m1,m2,x_vec,u_vec}); 
%matlabFunction(dx,"File","dx","Vars",[g,l1,l2,m1,m2,q1,q2,x0,y0,q1_d,q2_d,x0_d,y0_d,tau1,tau2,Fx,Fy]); 
clear dx; % Clear the variable dx to avoid conflict during simulation

% Create function KE.m for Kinetic Energy calculation
matlabFunction(K,"File","KE","Vars",{l1,lc1,l2,lc2,m1,m2,x_vec});
clear K; % Clear the variable K to avoid conflict during simulation

% Create function PE.m for Potential Energy calculation
matlabFunction(P,"File","PE","Vars",{g,l1,lc1,l2,lc2,m1,m2,x_vec});
clear P; % Clear the variable P to avoid conflict during simulation

% Create function A_matrix.m for A matrix of the linearized system
matlabFunction(A_mat,"File","A_matrix","Vars",{g,l1,lc1,l2,lc2,m1,m2,x_vec,u_vec});
clear A_mat;

% Create function B_matrix.m for B matrix of the linearized system
matlabFunction(B_mat,"File","B_matrix","Vars",{g,l1,lc1,l2,lc2,m1,m2,x_vec,u_vec});
clear B_mat;

