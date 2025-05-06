function dx = cubli_controlled(~, x, K, params)
    % Control input (LQR feedback)
    u = -K * x;

    % Extract states
    theta = x(1);
    omega = x(2);
    omega_w = x(3);

    % System parameters
    mh = params.mh;
    mw = params.mw;
    l = params.l;
    g = params.g;
    Jh = params.Jh;
    Jw = params.Jw;
    Km = params.Km;
    Cw = params.Cw;

    % Equations of motion (controlled)
    M = Jh + mw * l^2;
    torque_gravity = -mh * g * l * sin(theta);
    torque_control = Km * u; % Control input

    omega_dot = (torque_gravity + torque_control) / M;
    omega_w_dot = -omega_dot; % Reaction from wheel

    % State derivatives
    dx = [omega; omega_dot; omega_w_dot];
end
