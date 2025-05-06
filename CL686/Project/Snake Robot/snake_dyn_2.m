function Zdot = snake_dyn_2(t,Z,params,tspan,Ufk)
    % Z = [Phi,theta,px,py,Phi_dot,theta_dot,px_dot,py_dot]^T
    
    U = interp1(tspan,Ufk,t).';

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
