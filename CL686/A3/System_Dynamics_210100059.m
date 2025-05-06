 function Xdot = System_Dynamics_210100059(X, params, U, Fd)
    % Input  :  X - [h1,h2,h3,h4]^T, params are the physical parameters
    %           U - [v1,v2]^T, Fd - disturbance
    % Output :  Xdot - ODE equations for this system
    A1 = params.A1; A2 = params.A2; A3 = params.A3; A4 = params.A4;
    a1 = params.a1; a2 = params.a2; a3 = params.a3; a4 = params.a4;
    g = params.g;
    k1 = params.k1; k2 = params.k2;
    gamma1 = params.gamma1; gamma2 = params.gamma2; gamma3 = params.gamma3;

    Xdot(1) = -a1/A1*sqrt(2*g*X(1)) + a3/A1*sqrt(2*g*X(3)) + gamma1*k1/A1*U(1);
    Xdot(2) = -a2/A2*sqrt(2*g*X(2)) + a4/A2*sqrt(2*g*X(4)) + gamma2*k2/A2*U(2);
    Xdot(3) = -a3/A3*sqrt(2*g*X(3)) + (1-gamma2)*k2/A3*U(2) + gamma3/A3*Fd;
    Xdot(4) = -a4/A4*sqrt(2*g*X(4)) + (1-gamma1)*k1/A4*U(1) + (1-gamma3)/A4*Fd;
    Xdot = Xdot'; % because ODE45 apparently needs the eqns as a column vector
end