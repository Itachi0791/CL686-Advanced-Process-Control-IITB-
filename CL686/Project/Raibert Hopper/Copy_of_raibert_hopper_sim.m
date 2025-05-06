clc; clear; close all;

%% Parameters
params.g = 9.81; params.l0 = 2; params.k = 500; params.m = 1;
params.U = -pi/60; % Leg placement angle
dt = 0.001; T_fin = 5;
tspan = 0:dt:T_fin;
z0 = [0, 5, 3, 0]; % Initial state [x, xdot, y, ydot]

%% Initialize Simulation
% X_current = X0;
% t_current = 0;
% X_history = []; % Store all states
% t_history = []; % Store time history
% origin_history = []; % Store reference origin history
% phase_history = strings([]); % Store phase history
% phase = "flight"; % Start in flight phase

steps = 1;
fps = 30;

[z,t] = onestep_sim(z0,params,steps,fps);
plot(z(:,1),z(:,3))

%% Functions

function zdot = flight(~, z, params)
    % z = [x, xdot, y, ydot]^T
    g = params.g;
    zdot = [z(2); z(4); 0; -g];
end

function zdot = stance(~, z, params)
    x = z(1); y = z(3); %x & y position of com wrt ground
    l = sqrt(x^2+y^2);
    l0 = params.l0; k = params.k; m = params.m; g = params.g;
    F_spring = k*(l0-l);
    Fx_spring =  F_spring*(x/l);
    Fy_spring = F_spring*(y/l);
    Fy_gravity = m*g;
    xddot = (1/m)*(Fx_spring);
    yddot = (1/m)*(-Fy_gravity+Fy_spring);
    zdot = [z(2) xddot z(4) yddot]';
end

function [position, isterminal, direction] = touchdown(~, z, params)
    l0 = params.l0; theta = params.U;
    position = z(3) - l0 * cos(theta);
    isterminal = 1;
    direction = -1;
end

function [position, isterminal, direction] = liftoff(~, z, params)
    l = sqrt(z(1)^2+z(3)^2);
    l0 = params.l0;
    position = l - l0;
    isterminal = 1;
    direction = 1;
end

function [gstop, isterminal,direction]=apex(~,z,~)
    gstop = z(4) - 0; %ydot is 0;
    direction = 0; %negative direction goes from + to -
    isterminal = 1;  %1 is stop the integration
end

function animate(t_all,z_all,params,fps)
    l0 = params.l0;
   
    %%% interpolate for animation %%
    [t_interp,z_interp] = loco_interpolate(t_all,z_all,fps);
    
    %%%%% prepare for animation %%%%%%%
    [mm,~] = size(z_interp);
    min_xh = min(z_interp(:,1)); max_xh = max(z_interp(:,1)); 
    dist_travelled = max_xh - min_xh;
    camera_rate = dist_travelled/mm;
    
    window_xmin = -1.0*l0; window_xmax = l0;
    window_ymin = -0.1; window_ymax = 1.9*l0;
    
    axis('equal')
    axis([window_xmin window_xmax window_ymin window_ymax])
    axis off
    set(gcf,'Color',[1,1,1])
    
    
    figure(1);
    for i=1:length(t_interp)
       
        plot(z_interp(i,1),z_interp(i,3),'ro','MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',20); %com
        line([-1 max(z_interp(:,1))+1],[0 0],'Linewidth',2,'Color','black'); %ground
        line([z_interp(i,1) z_interp(i,5)],[z_interp(i,3) z_interp(i,6)],'Linewidth',4,'Color',[0 0.8 0]); %leg
         
        window_xmin = window_xmin + camera_rate;
        window_xmax = window_xmax + camera_rate;
        axis('equal')
        axis off
        axis([window_xmin window_xmax window_ymin window_ymax])
    
        pause(0.05);
    end
end


function [t_interp,z_interp] = loco_interpolate(t_all,z_all,fps)
    [~,n] = size(z_all);
    t_interp = linspace(t_all(1),t_all(end),fps*(t_all(end)-t_all(1)));
    
    for i=1:n
        z_interp(:,i) = interp1(t_all,z_all(:,i),t_interp);
    end
    t_interp = t_interp';
end

function [z,t]= onestep_sim(z0,params,steps,fps)  
    l0 = params.l0; theta = params.U;
    flag = 1;
    if nargin<2
        error('need more inputs to onestep');
    elseif nargin<3
        flag = 0; %send only last state, for root finder and jacobian
        steps = 1;
        fps = 50;
    end
    
    %x0 = 0; x0dot = z0(1);  
    %y0 = z0(2); y0dot = 0;
    
    %z0 = [x0 x0dot y0 y0dot];
    
    t0 = 0; 
    dt = 5; %might need to be changed based on time taken for one step
    %time_stamps = 100;
    t_ode = t0;
    z_ode = [z0 ...
             z0(1)+l0*sin(theta) ...
             z0(3)-l0*cos(theta)];
    
    for i=1:steps
        
        %%% apex to ground %%%
        options1 = odeset('Abstol',1e-13,'Reltol',1e-13,'Events',@touchdown);
        tspan = linspace(t0,t0+dt,dt*1000);
        [t_temp1,z_temp1]=ode113(@flight,tspan,z0,options1,params);
        [t_temp1,z_temp1] = loco_interpolate(t_temp1,z_temp1,10*fps);
    
        z_temp1 = [z_temp1 ...
                   z_temp1(1:end,1)+l0*sin(theta) ...
                   z_temp1(1:end,3)-l0*cos(theta) ...
                  ];
        t0 = t_temp1(end);
        z0(1:4) = z_temp1(end,1:4);
        x_com = z0(1); %save the x position for future
        z0(1) = -l0*sin(theta); %relative distance wrt contact point because of non-holonomic nature of the system
        x_foot = x_com + l0*sin(theta); 
        y_foot = 0;
       
        %%% stance phase %%%
        tspan = linspace(t0,t0+dt,dt*2000);
        options2 = odeset('Abstol',1e-13,'Reltol',1e-13,'Events',@liftoff);
        [t_temp2,z_temp2]=ode113(@stance,tspan,z0,options2,params);
        [t_temp2,z_temp2] = loco_interpolate(t_temp2,z_temp2,10*fps);
        
        z_temp2(:,1) = z_temp2(:,1) + x_com + l0*sin(theta); %absolute x co-ordinate
        z_temp2 = [z_temp2, ...
              x_foot*ones(length(z_temp2),1) y_foot*ones(length(z_temp2),1)]; %the distal end of leg is 0 when touching the ground.
        t0 = t_temp2(end);
        z0(1:4) = z_temp2(end,1:4);
        
        %%% ground to apex
        tspan = linspace(t0,t0+dt,dt*1000);
        options3 = odeset('Abstol',1e-13,'Reltol',1e-13,'Events',@apex);
        [t_temp3,z_temp3]=ode113(@flight,tspan,z0,options3,params);
        [t_temp3,z_temp3] = loco_interpolate(t_temp3,z_temp3,10*fps);
    
        
         z_temp3 = [z_temp3 ...
                   z_temp3(1:end,1)+l0*sin(theta) ...
                   z_temp3(1:end,3)-l0*cos(theta) ...
                  ];
        t0 = t_temp3(end);
        z0(1:4) = z_temp3(end,1:4);
        
        %%%%% Ignore time stamps for heelstrike and first integration point
        t_ode = [t_ode; t_temp1(2:end); t_temp2(2:end);  t_temp3(2:end)];
        z_ode = [z_ode; z_temp1(2:end,:); z_temp2(2:end,:); z_temp3(2:end,:)];
        
    end
    
    z = [z0(2) z0(3)];
    
    if flag==1
       z=z_ode;
       t=t_ode;
    end
end
