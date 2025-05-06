% clc; clear; close all;
% 
% %% Parameters
% params.g = 9.81; params.l0 = 3; params.k = 1000; params.m = 1;
% U = 0; % Leg placement angle
% dt = 0.02;T_fin = 10;
% tspan = 0:dt:T_fin; 
% X0 = [4, 5, 0, 10]; % Initial state [x, z, xdot, zdot]
% 
% %% Flight Dynamics
% function Xdot = flight_dynamics(~, X, params)
%     % X = [x, z, xdot, zdot]^T
%     g = params.g;
%     Xdot = [X(3); X(4); 0; -g];
% end
% 
% %% Stance Dynamics
% % [t,X_rot] = ode45(@(t,X)stance_dynamics(t, X, params),0:0.02:1,X_stance0);
% % X_rot = X_rot';
% % X_cart = [-X_rot(1,:).*sin(X_rot(2,:));...
% %             X_rot(1,:).*cos(X_rot(2,:));...
% %            -X_rot(1,:).*cos(X_rot(2,:)).*X_rot(4,:)-X_rot(3,:).*sin(X_rot(2,:));
% %             X_rot(3,:).*cos(X_rot(2,:))-X_rot(1,:).*sin(X_rot(2,:)).*X_rot(4,:)];
% % figure;
% % plot(X_cart(1,:),X_cart(2,:),'linewidth',2)
% % axis equal;axis([-2 2 0 2])
% function Xdot = stance_dynamics(~, X, params)
%     % X = [r,theta,rdot,thetadot]^T
%     g = params.g; l0 = params.l0; m = params.m; k = params.k;
%     Xdot = [X(3); X(4);
%         X(1) * X(4).^2 - g * cos(X(2)) + k/m * (l0 - X(1));
%         g / X(1) * sin(X(2)) + 2 * X(3) / X(1) * X(4)];
% end
% 
% %% Touchdown Event
% function [position, isterminal, direction] = touchdown(~, X, params, U)
%     l0 = params.l0;
%     position = X(2) - l0 * cos(U); % Detect when the foot reaches ground
%     isterminal = 1; % Stop integration
%     direction = -1; % Detect downward crossing
% end
% 
% %% Liftoff Condition
% function [position, isterminal, direction] = liftoff(~, X, params)
%      l0 = params.l0;
%      position = X(1) - l0; % Lift-off when leg reaches neutral length
%      isterminal = 1; % Stop integration
%      direction = 1; % Detect downward crossing
% end
% 
% %% Simulate Flight Phase
% 
% options1 = odeset("AbsTol",1e-10,"RelTol",1e-10, ...
%                  "Events", @(t, X) touchdown(t, X, params, U));
% options2 = odeset("AbsTol",1e-10,"RelTol",1e-10, ...
%                  "Events", @(t, X) liftoff(t, X, params));
% 
% [t_flight, X_flight] = ode45(@(t, X) flight_dynamics(t, X, params), tspan, X0, options1);
% t_flight = t_flight';X_flight = X_flight';
% new_origin = [X_flight(1,end) + params.l0 * sin(U);X_flight(2,end) - params.l0 * cos(U);0;0];
% X_stance0 = fl_to_st(X_flight(:,end),U,params);
% [t_stance, X_stance] = ode45(@(t, X) stance_dynamics(t, X, params), t_flight(end)+tspan, X_stance0, options2);
% t_stance = t_stance';X_stance = X_stance';
% X_stance_cart = [new_origin(1)-X_stance(1,:).*sin(X_stance(2,:));...
%             X_stance(1,:).*cos(X_stance(2,:));...
%            -X_stance(1,:).*cos(X_stance(2,:)).*X_stance(4,:)-X_stance(3,:).*sin(X_stance(2,:));
%             X_stance(3,:).*cos(X_stance(2,:))-X_stance(1,:).*sin(X_stance(2,:)).*X_stance(4,:)];
% 
% 
% X_flight_1 = st_to_fl(X_stance(:,end),new_origin);
% [t_flight_2, X_flight_2] = ode45(@(t, X) flight_dynamics(t, X, params), t_stance(end)+t_flight(end)+tspan, X_flight_1, options1);
% t_flight_2 = t_flight_2';X_flight_2 = X_flight_2';
% 
% X = [X_flight,X_stance_cart,X_flight_2]; t = [t_flight,t_stance,t_flight_2];
% %% Animation
% figure;
% hold on; %axis equal;
% xlim([-1, 10]); ylim([0, 11]);
% hopper_mass = plot(X0(1), X0(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
% hopper_leg = plot([X0(1), X0(1)], [0, X0(2)], 'b-', 'LineWidth', 2);
% 
% for i = 1:length(t_flight)
%     % Update hopper position
%     set(hopper_mass, 'XData', X_flight(1,i), 'YData', X_flight(2,i));
%     set(hopper_leg, 'XData', [X_flight(1,i), X_flight(1,i) + params.l0 * sin(U)], ...
%                     'YData', [X_flight(2,i), X_flight(2,i) - params.l0 * cos(U)]);
%     pause(0.02);
% end
% for i = 1:length(t_stance)
%     % Update hopper position
%     set(hopper_mass, 'XData', X_stance_cart(1,i), 'YData', X_stance_cart(2,i));
%     set(hopper_leg, 'XData', [X_stance_cart(1,i), X_stance_cart(1,i) + X_stance(1,i).*sin(X_stance(2,i))], ...
%                     'YData', [X_stance_cart(2,i), X_stance_cart(2,i) - X_stance(1,i).* cos(X_stance(2,i))]);
%     pause(0.02);
% end
% 
% for i = 1:length(t_flight_2)
%     % Update hopper position
%    set(hopper_mass, 'XData', X_flight_2(1,i), 'YData', X_flight_2(2,i));
%     set(hopper_leg, 'XData', [X_flight_2(1,i), X_flight_2(1,i) + params.l0 * sin(U)], ...
%                     'YData', [X_flight_2(2,i), X_flight_2(2,i) - params.l0 * cos(U)]);
% 
%     pause(0.02);
% end
% 
% function X_stance = fl_to_st(X_flight,U,params)
%     rdot = X_flight(4).*cos(U) - X_flight(3).*sin(U);
%     thetadot = -1/params.l0*(X_flight(3).*cos(U) + X_flight(4).*sin(U));
%     X_stance = [params.l0;U;rdot;thetadot];
% end
% 
% function X_flight = st_to_fl(X_stance,origin)
%     r = X_stance(1);
%     theta = X_stance(2);
%     r_dot = X_stance(3);
%     theta_dot = X_stance(4);
%     xdot = -r_dot.*sin(theta) - r*theta_dot*cos(theta);
%     ydot = r_dot.*cos(theta) - r*theta_dot*sin(theta);
%     X_flight = [origin(1) - r * sin(theta);origin(2) + r * cos(theta);xdot;ydot];
% end
% 
% 

clc; clear; close all;

%% Parameters

params.g = 10; params.l0 = 1; params.k = 100; params.m = 1;
U = -10*(pi/180); % Leg placement angle
U_history = [];
dt = 0.001; T_fin = 5;
tspan = 0:dt:T_fin;
X0 = [0, 1.182441943200824, 1.017450873821607, 0]; % Initial state [x, z, xdot, zdot]

  
options1 = odeset("AbsTol",1e-10,"RelTol",1e-10, ...
                 "Events", @(t, X) touchdown(t, X, params, U));
options2 = odeset("AbsTol",1e-10,"RelTol",1e-10, ...
                 "Events", @(t, X) liftoff(t, X, params));

%% Initialize Simulation
X_current = X0;
t_current = 0;
X_history = []; % Store all states
t_history = []; % Store time history
origin_history = []; % Store reference origin history
phase_history = strings([]); % Store phase history
phase = "flight"; % Start in flight phase

%% Simulation Loop
while t_current < T_fin
   
    U = U;

    if phase == "flight"
        % Simulate Flight Phase
        [t_flight, X_flight,te,xe,ie] = ode45(@(t, X) flight_dynamics(t, X, params), tspan, X_current, options1);
        t_flight = t_flight'; X_flight = X_flight';
        disp(ie)
        % Store results
        X_history = [X_history, X_flight];
        t_history = [t_history, t_current + t_flight];
        phase_history = [phase_history,repmat("flight",1,size(X_flight,2))];
        origin_history = [origin_history,repmat([0;0],1,size(X_flight,2))];
        U_history = [U_history, repmat(U, 1, size(X_flight, 2))];

        if isempty(t_flight) % If no touchdown, stop simulation
            break;
        end
       
        % Transition to Stance Phase
        new_origin = [X_flight(1,end) + params.l0 * sin(U); X_flight(2,end) - params.l0 * cos(U)];
        X_stance0 = fl_to_st(X_flight(:,end), U, params);
        X_current = X_stance0;
        phase = "stance";
        t_current = t_current + t_flight(end);

    elseif phase == "stance"
        % Simulate Stance Phase
        [t_stance, X_stance] = ode45(@(t, X) stance_dynamics(t, X, params), tspan, X_current, options2);
        t_stance = t_stance'; X_stance = X_stance';

        % Convert to Cartesian Coordinates
        X_stance_cart = [new_origin(1) - X_stance(1,:) .* sin(X_stance(2,:)); ...
                         X_stance(1,:) .* cos(X_stance(2,:)); ...
                        -X_stance(1,:) .* cos(X_stance(2,:)) .* X_stance(4,:) - X_stance(3,:) .* sin(X_stance(2,:)); ...
                         X_stance(3,:) .* cos(X_stance(2,:)) - X_stance(1,:) .* sin(X_stance(2,:)) .* X_stance(4,:)];

        

        % Store results
        X_history = [X_history, X_stance_cart];
        t_history = [t_history, t_current + t_stance];
        phase_history = [phase_history,repmat("stance",1,size(X_stance,2))];
        origin_history = [origin_history,repmat(new_origin,1,size(X_stance,2))];
        U_history = [U_history, repmat(U, 1, size(X_stance, 2))];

        U = X_stance(2,end);
        if isempty(t_stance) % If no liftoff, stop simulation
            break;
        end

        % Transition to Flight Phase
        X_current = st_to_fl(X_stance(:,end), new_origin);
        phase = "flight";
        t_current = t_current + t_stance(end);

       

    end


end

%% Animation
figure;
hold on;
plot(-50:50,zeros(101,1),'k')
xlim([-20, 50]); ylim([-5, 11]);
hopper_mass = plot(X0(1), X0(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
hopper_leg = plot([X0(1), X0(1)], [0, X0(2)], 'b-', 'LineWidth', 2);

for i = 1:length(t_history)
    % Update hopper position
    set(hopper_mass, 'XData', X_history(1,i), 'YData', X_history(2,i));
    
    % Determine phase (flight or stance)
    if phase_history(i) == "flight"  % Flight phase
        foot_x = X_history(1,i) + params.l0 * sin(U_history(i));
        foot_y = X_history(2,i) - params.l0 * cos(U_history(i));
    else  % Stance phase
        foot_x = origin_history(1,i);
        foot_y = origin_history(2,i);
    end
    foot_history(:,i) = [foot_x;foot_y];

    % Update leg position
    set(hopper_leg, 'XData', [X_history(1,i), foot_x], ...
                    'YData', [X_history(2,i), foot_y]);
    if i == length(t_history)
        pause(dt);
    end
    pause(0.001);%pause(t_history(i+1)-t_history(i));
end

%% Functions

function Xdot = flight_dynamics(~, X, params)
    % X = [x, z, xdot, zdot]^T
    g = params.g;
    Xdot = [X(3); X(4); 0; -g];
end

function Xdot = stance_dynamics(~, X, params)
    % X = [r,theta,rdot,thetadot]^T
    g = params.g; l0 = params.l0; m = params.m; k = params.k;
    Xdot = [X(3); X(4);
        X(1) * X(4).^2 - g * cos(X(2)) + k/m * (l0 - X(1));
        g / X(1) * sin(X(2)) + 2 * X(3) / X(1) * X(4)];
end

function [position, isterminal, direction] = touchdown(~, X, params, U)
    l0 = params.l0;
    position = X(2) - l0 * cos(U);
    isterminal = 1;
    direction = -1;
end

function [position, isterminal, direction] = liftoff(~, X, params)
    l0 = params.l0;
    position = X(1) - l0;
    isterminal = 1;
    direction = 1;
end

function X_stance = fl_to_st(X_flight, U, params)
    rdot = X_flight(4) .* cos(U) - X_flight(3) .* sin(U);
    thetadot = -1/params.l0 * (X_flight(3) .* cos(U) + X_flight(4) .* sin(U));
    X_stance = [params.l0; U; rdot; thetadot];
end

function X_flight = st_to_fl(X_stance, origin)
    r = X_stance(1);
    theta = X_stance(2);
    r_dot = X_stance(3);
    theta_dot = X_stance(4);
    xdot = -r_dot .* sin(theta) - r * theta_dot * cos(theta);
    ydot = r_dot .* cos(theta) - r * theta_dot * sin(theta);
    X_flight = [origin(1) - r * sin(theta); origin(2) + r * cos(theta); xdot; ydot];
end


