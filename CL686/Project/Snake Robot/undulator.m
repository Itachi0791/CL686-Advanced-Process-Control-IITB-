function U = undulator(t,Z,params,Kp,Kd,Phi_dd_ref_fun, Phi_d_ref_fun, Phi_ref_fun)
    U = Phi_dd_ref_fun(t) + Kp * (Phi_ref_fun(t) - Z(1:params.N_l-1,1)) + Kd * (Phi_d_ref_fun(t) - Z(params.N_l+3: 2*params.N_l+1,1));
end