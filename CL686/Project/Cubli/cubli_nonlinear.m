function dx = cubli_nonlinear(~, x, params)
    % Extract states
    theta = x(1); % Tilt angle
    omega = x(2); % Angular velocity
    omega_w = x(3); % Wheel velocity

    % System parameters
    mh = params.mh;
    mw = params.mw;
    l = params.l;
    g = params.g;
    Jh = params.Jh;
    Jw = params.Jw;
    Km = params.Km;
    Cw = params.Cw;

    % Equations of motion (nonlinear)
    M = Jh + mw * l^2;
    torque_gravity = -mh * g * l * sin(theta);
    torque_wheel = -Km * omega_w - Cw * omega_w;

    omega_dot = (torque_gravity + torque_wheel) / M;
    omega_w_dot = -omega_dot; % Reaction from wheel

    % State derivatives
    dx = [omega; omega_dot; omega_w_dot];
end
