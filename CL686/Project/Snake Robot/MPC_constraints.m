function [c,ceq] = MPC_constraints(Ts,Np,Z0,params,Ufk)

     % Z = [Phi,theta,px,py,Phi_dot,theta_dot,px_dot,py_dot]^T
    N_l = params.N_l;
    Z = zeros(Np,length(Z0)); 
    Z(1,:) = Z0; 

    for k = 1:Np-1
        Z(k+1,:) = Z(k,:) + Ts*snake_dyn(Z(k,:)',params,Ufk(k,:)')';
    end

    c = zeros((N_l-1)*Np,1);

    for i=1:Np
        for j=1:N_l-1
            c(j+(i-1)*(N_l-1)) = (Z(i,j)-0.05)*(Z(i,j)+0.05);
        end
    end

    ceq = [];
end