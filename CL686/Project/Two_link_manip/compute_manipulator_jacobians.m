function compute_manipulator_jacobians(m1,m2,l1,l2,g)
    % Define symbolic variables for states and inputs
    syms theta1 theta2 theta1_dot theta2_dot u1 u2 real
    
    I1 = m1*l1^2/3;     % Inertia of link 1
    I2 = m2*l2^2/3;     % Inertia of link 2
    
    % Inertia matrix H(q)
    H = [I1 + m1*l1^2/4 + I2 + m2*(l2^2/4 + l1^2 + l1*l2*cos(theta2)), ...
          I2 + m2*(l2^2/4 + l1*l2*cos(theta2)/2); ...
          I2 + m2*(l2^2/4 + l1*l2*cos(theta2)/2), ...
          I2 + m2*l2^2/4];
    
    % Coriolis matrix C(q, q_dot)
    C = [-m2*l1*l2*sin(theta2)*theta2_dot, -m2*l1*l2*sin(theta2)*theta2_dot/2; ...
          m2*l1*l2*sin(theta2)*theta1_dot/2, 0];
    
    % Gravity vector G(q)
    G = [m1*g*l1*cos(theta1)/2 + m2*g*(l1*cos(theta1) + l2/2*cos(theta1 + theta2)); ...
          m2*g*l2/2*cos(theta1 + theta2)];
    
    % Joint accelerations q_ddot = H^{-1}(u - C*q_dot - G)
    u = [u1; u2];
    q_dot = [theta1_dot; theta2_dot];
    q_ddot = H \ (u - C*q_dot - G);
    
    % State-space dynamics: x_dot = f(x, u)
    x_dot = [theta1_dot; theta2_dot; q_ddot];
    
    % State vector: x = [theta1; theta2; theta1_dot; theta2_dot]
    x = [theta1; theta2; theta1_dot; theta2_dot];

    % Output dynamics
    y = [l1*cos(theta1)+l2*cos(theta2);l1*sin(theta1)+l2*sin(theta2)];
    
    % Compute Jacobians A = df/dx, B = df/du
    A_matrix = jacobian(x_dot, x);
    B_matrix = jacobian(x_dot, u);
    C_matrix = jacobian(y,x);
    matlabFunction(A_matrix,"File","A_matrix", 'Vars', {[theta1; theta2; theta1_dot; theta2_dot], [u1; u2]});
    matlabFunction(B_matrix,"File","B_matrix", 'Vars', {[theta1; theta2; theta1_dot; theta2_dot], [u1; u2]});
    matlabFunction(C_matrix,"File","C_matrix", 'Vars', {[theta1; theta2; theta1_dot; theta2_dot], [u1; u2]});
end

