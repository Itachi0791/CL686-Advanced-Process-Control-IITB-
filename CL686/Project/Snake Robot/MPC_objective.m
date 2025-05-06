function obj_val = MPC_objective(Ts,Np,Z0,params,Ufk)

     % Z = [Phi,theta,px,py,Phi_dot,theta_dot,px_dot,py_dot]^T
    t_dur = 0:Ts:Ts*(Np-1);
    N = params.N_l;

    Z = zeros(length(t_dur),length(Z0)); 
    Z(1,:) = Z0; 
    obj_val = 0;

    for k = 1:length(t_dur)-1
        Z(k+1,:) = Z(k,:) + Ts*snake_dyn(Z(k,:)',params,Ufk(k,:)')';
        theta= Z(k,N);
        px_dot = Z(k,2*N+3);
        py_dot = Z(k,2*N+4);
        v_t = cos(theta)*px_dot + sin(theta)*py_dot;
        obj_val = obj_val-v_t + 0.05 * Z(k,1:N-1)*Z(k,1:N-1)'; % Objective Function Value
    end

    % 
    % theta_end = Z(end,N);
    % px_dot_end = Z(end,2*N+3);
    % py_dot_end = Z(end,2*N+4);
    % v_t_end = cos(theta_end)*px_dot_end + sin(theta_end)*py_dot_end;
    % 
    % obj_val = -v_t_end; % Objective Function Value
end




