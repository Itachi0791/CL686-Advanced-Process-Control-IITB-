function [K, A, B] = cubli_lqr(params)
    % Extract parameters
    mh = params.mh;
    mw = params.mw;
    l = params.l;
    g = params.g;
    Jh = params.Jh;
    Jw = params.Jw;
    Km = params.Km;
    Cw = params.Cw;

    % Linearized System Matrices
    M = Jh + mw * l^2;
    A = [0 1 0; mh * g * l / M 0 -1 / M; 0 0 0];
    B = [0; Km / M; -Km / Jw];

    % LQR Control Design
    Q = diag([10, 1, 0.1]); % Penalize angle and velocity
    R = 1; % Control effort penalty
    K = lqr(A, B, Q, R);
end
