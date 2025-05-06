%% Loco-manipulator simulation in MATLAB

clc;clear;
close all;

% Link 1 = Mobile Base
% Link 2 = Manipulator Link 1
% Link 3 = Manipulator Link 2

%% System Parameters

m1 = 1; % Mass of link 1
m2 = 1; % Mass of link 2
m3 = 1; % Mass of link 3
l1 = 1; % Length of link 1
l2 = 1; % Length of link 2
l3 = 1; % Length of link 3
g = 0; % Acceleration due to gravity 
I1 = (m1*l1^2)/12; % Moment of inertia of link 1 about COM 
I2 = (m2*l2^2)/12; % Moment of inertia of link 2 about COM
I3 = (m3*l3^2)/12; % Moment of inertia of link 3 about COM
lc1 = (l1/2); % Location of COM of link 1
lc2 = (l2/2); % Location of COM of link 2
lc3 = (l3/2); % Location of COM of link 3

dt = 0.01; % Time step
Tmax = 30; % Maximum time

tspan = 0:dt:Tmax; % Time duration for the simulation

x0=[0;0;0;0;0;0;0;0;0;0;0]; % Initial state

%% Actuator Inputs on the System

u1 = @(t,x)10*(pi/3 - x(1)) + 2*(-x(6)); % Torque on mobile base
u2 = @(t,x)10*(pi/2 - x(2)) + 2*(-x(7)); % Torque on manipulator link 1 
u3 = @(t,x)10*(pi/2 - x(3)) + 2*(-x(8)); % Torque on manipulator link 2 
Ft = @(t,x)0; % Thrust on the mobile base

%% Numerical DE solving using ode45

options=odeset('abstol',1e-9,'reltol',1e-9);
[t,x] = ode45(@(t,x)dx(t,g,l1,l2,l3,m1,m2,m3,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),u1(t,x),u2(t,x),u3(t,x),Ft(t,x)),tspan,x0,options); 

%% Calculation of different variables (acceleration, energy, torque)

% Vector storing the kinetic energy at different time instances
K = zeros(length(x),1); 
% Vector storing the potential energy at different time instances
P = zeros(length(x),1); 
% Vector storing the total energy at different time instances
T = zeros(length(x),1); 
% Vector storing the actuator work at different time instances
W_ac = zeros(length(x),1); 

for j=1:length(x)
    % Kinetic energy
    K(j,1) = KE(l1,l2,l3,m1,m2,m3,x(j,1),x(j,2),x(j,3),x(j,4),x(j,5),x(j,6),x(j,7),x(j,8),x(j,9),x(j,10));
    % Potential energy
    P(j,1) = PE(g,l1,l2,l3,m1,m2,m3,x(j,1),x(j,2),x(j,3),x(j,4),x(j,5),x(j,6),x(j,7),x(j,8),x(j,9),x(j,10));
    % Total Energy
    T(j,1) = K(j,1) + P(j,1);
    % Actuator Work
    W_ac(j,1) = x(j,end); 
end

%% Animation

set(gcf,'defaultAxesTickLabelInterpreter','latex'); 

for i=1:10:length(t)
    x1 = x(i,4) + l1 * cos(x(i,1));
    y1 = x(i,5) + l1 * sin(x(i,1));
    x2 = x1 + l2 * cos(x(i,1) + x(i,2));
    y2 = y1 + l2 * sin(x(i,1) + x(i,2));
    x3 = x2 + l3 * cos(x(i,1) + x(i,2) + x(i,3));
    y3 = y2 + l3 * sin(x(i,1) + x(i,2) + x(i,3));

    scatter (x(i,4),x(i,5), 200, "black", "filled")
    hold on
    plot([x(i,4),x1],[x(i,5),y1],'r-o','LineWidth',5)
    hold on
    plot([x1,x2],[y1,y2],'b-o','LineWidth',5)
    hold on
    plot([x2,x3],[y2,y3],'magenta','LineWidth',5)

    grid on, set(gca, 'FontSize',15)
    axis([-10,10,-10,10])
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    title('Animation','Interpreter','latex')
    time = Tmax*(i-1)/(length(t)-1);
    subtitle("Time = " + time + " sec",'Interpreter','latex')
    hold off

    pause(0.01)

end

%% Energy Check

set(gcf,'defaultAxesTickLabelInterpreter','latex'); 

plot(t, T - W_ac - T(1,1),"LineWidth",1.5,"Color","Blue") % Kinetic energy vs time
xlabel("$t$[s]","FontSize",15,"Interpreter","latex")
ylabel("Energy Difference [J]","FontSize",15,"Interpreter","latex")
title("Energy Check","FontSize",15,"Interpreter","latex")
set(gca,'FontSize',15,'TickLabelInterpreter','latex')
