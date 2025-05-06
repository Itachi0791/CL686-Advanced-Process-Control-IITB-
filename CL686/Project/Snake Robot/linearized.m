%% Deriving Matrices for Linearized System
function linearized(N_l)

    A = zeros(N_l-1,N_l); D = A;
    for i = 1:N_l-1
        A(i,i:i+1) = [1 1];
        D(i,i:i+1) = [1 -1];
    end
    e_bar = ones(N_l-1,1);

    Phi = sym("Phi",[N_l - 1,1]);
    Phi_dot = sym("Phi_dot",[N_l - 1,1]);
    U = sym("U",[N_l - 1, 1]);
    syms theta theta_dot p_x p_y px_dot py_dot 
    syms m cp cn ct lambda1 lambda2

    Z = [Phi; theta; p_x; p_y; Phi_dot; theta_dot; px_dot; py_dot];

    v_t = cos(theta)*px_dot + sin(theta)*py_dot;
    v_n = -sin(theta)*px_dot + cos(theta)*py_dot;
    Phi_ddot = (-cn/m)*Phi_dot + (cp/m)*v_t*A*D'*Phi+(1/m)*(D*D')*U;
    theta_ddot = (-lambda1)*theta_dot + (lambda2/(N_l-1))*v_t*e_bar'*Phi;
    v_t_dot = (-ct/m)*v_t + (2*cp/(N_l*m))*v_n*e_bar'*Phi - (cp/(N_l*m))*Phi'*A*pinv(D)*Phi_dot;
    v_n_dot = (-cn/m)*v_n + (2*cp/(N_l*m))*v_t*e_bar'*Phi;
    px_ddot = v_t_dot*cos(theta) -v_n_dot*sin(theta) -v_t*sin(theta)*theta_dot -v_n*cos(theta)*theta_dot;
    py_ddot = v_t_dot*sin(theta) +v_n_dot*cos(theta) +v_t*cos(theta)*theta_dot -v_n*sin(theta)*theta_dot;

    f = [Phi_dot;theta_dot;px_dot;py_dot;Phi_ddot;theta_ddot;px_ddot;py_ddot];

    A_matrix = jacobian(f,Z);
    B_matrix = jacobian(f,U);

    matlabFunction(A_matrix,"File","A_matrix","Vars",{Z,m,cp,cn,ct,lambda1,lambda2});
    matlabFunction(B_matrix,"File","B_matrix");
end