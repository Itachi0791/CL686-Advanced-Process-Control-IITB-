function [Sx,Su,Sdel,H,F] = MPC_matrices(Phi_mat,Gamma_mat,p,Wx,Wu,x_k,Del_k,Zs)
    n = size(Phi_mat,1); m = size(Gamma_mat,2);
    Sx = zeros(n*p,n); Su = zeros(n*p, m*p); Sdel = zeros(n*p,n*p); % Matrices in equality constraint, X = Sx*x(k) + Su*Uf,k + Sdelta*Deltak
    for i = 1:p
        Sx(n*(i-1)+1:n*i,1:n) = Phi_mat^i; 
        for j = 1:i
            Su(n*(i-1)+1:n*i,m*(i-j)+1:m*(i-j)+m) = Phi_mat^(j-1)*Gamma_mat;
            Sdel(n*(i-1)+1:n*i,n*(i-1)+1:n*i) = Phi_mat^(j-1);
        end
    end
    H = 2*(Su'*Wx*Su + Wu); H = (H + H')/2;
    F = 2*(Wx*Su)'*(Sx*x_k+Sdel*Del_k-Zs);
end